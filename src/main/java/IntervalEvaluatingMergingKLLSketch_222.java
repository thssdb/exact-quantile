import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class IntervalEvaluatingMergingKLLSketch_222 {
    int dataType = 1;
    static int pageN = 8192;
    static int N = 55000000/pageN*pageN, pageNum=N/pageN; // CHECK IT
    public static int TEST_CASE=32; // CHECK IT
    boolean TEST_FULL=true; // CHECK IT
    static double[] a;
    static KLLSketchForQuantile[] KLLArr;
//    String time_result="";
    static ArrayList<String> err_result = new ArrayList<>();
    static ArrayList<String> time_result = new ArrayList<>();
    boolean show_time = false,show_err=true;
    int RESULT_LINE=0;
    Random random = new Random(233);

    public void prepareA(int dataType) throws IOException {
        if(a==null)a = new double[N];
        this.dataType = dataType;

        if(dataType==0){
            for(int i=0;i<N;i++)a[i]=
//                Math.pow(-1,random.nextInt(2))*Math.pow(10.0,(2*Math.pow(random.nextDouble(),2)-1)*300);
                (double)(i/10000000);
//            System.out.println("\t\t???a????\t"+dataToLong(a[0])+"\t\t+"+longToResult(dataToLong(a[0])));
        }
        if(dataType==1){
            BufferedReader reader = new BufferedReader(new FileReader(new File("1_bitcoin.csv")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN= 0;
            while((line=reader.readLine())!=null){
                a[cntN++] = Double.parseDouble(line);
                if(cntN==N)return;
            }
        }
        if(dataType==2){
            BufferedReader reader = new BufferedReader(new FileReader(new File("2_SpacecraftThruster.txt")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN= 0;
            while((line=reader.readLine())!=null){
                a[cntN++] = Double.parseDouble(line);
                if(cntN==N)return;
            }
        }
        if(dataType==3){
            BufferedReader reader = new BufferedReader(new FileReader(new File("3_taxipredition8M.txt")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN= 0;
            while((line=reader.readLine())!=null){
                a[cntN++] = Double.parseDouble(line);
                if(cntN==N)return;
            }
        }
        if(dataType==4){
            BufferedReader reader = new BufferedReader(new FileReader(new File("4_wh.csv")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN= 0;
            while((line=reader.readLine())!=null){
                a[cntN++] = Double.parseDouble(line);
                if(cntN==N)return;
            }
        }
    }
    public void prepareWorker(int maxSeriesByte){
        KLLArr = new KLLSketchForQuantile[pageNum];
        int enoughMemByte = pageN*10;
        for(int i=0;i<pageNum;i++) {
            LongKLLSketch worker = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);
            for (int j = 0; j<pageN; j++) worker.update(dataToLong(a[i*pageN+j]));
            worker.compactBeforeSerialization();
            KLLArr[i] = worker;
//            if(i==0)worker.show();
        }
    }
    public int getValueActualRank(double[]sortedA,int queryN, double v){ // number of elements <= v
        int L=0,R=queryN-1;
        while(L<R){
            int mid=(L+R+1)>>>1;
            if(v<sortedA[mid])R=mid-1;
            else L=mid;
        }
        return L;
    }
    public int getValueLessThan(double[]sortedA,int queryN, double v){ // number of elements <= v
        int L=0,R=queryN-1;
        while(L<R){
            int mid=(L+R+1)>>>1;
            if(sortedA[mid]<v)L=mid;
            else R=mid-1;
        }
        return sortedA[L]<v?L:L-1;
    }
    public int getDeltaRank(double[]sortedA,int queryN, double v,int targetRank){
        int rank_L = getValueLessThan(sortedA,queryN,v)+1;
        int rank_R = getValueActualRank(sortedA,queryN,v);
//        System.out.println("\t\t\t"+targetRank+"\t\tresultLR:"+rank_L+"..."+rank_R+"\t\tresV:"+v);
        if(targetRank>=rank_L&&targetRank<=rank_R)return 0;
        else return targetRank<rank_L?(targetRank-rank_L):(targetRank-rank_R);
    }



    private long dataToLong(double data)  {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }

    private double longToResult(long result)  {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }

    public void testMergeError(int queryN, int maxMemoryByte, int maxSeriesByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time=0, merge_time=0;
        double err_full=0,err_merge=0;
        double lv_full=0,lv_merge=0;
        double MMP_merge = 0,MMP_full = 0;
        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE*10];
        int[] RR = new int[TEST_CASE*10];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE*10; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }

        for(int T=0,cntTestId=0;T<TEST_CASE;T++,cntTestId++){
            if((cntTestId&1)==0)//CHECK IT!
                prepareWorker(maxSeriesByte);
            int L = LL[cntTestId], R = RR[cntTestId];
            int pageL = (L+pageN-1)/pageN, pageR = R/pageN;
            int posL = pageL*pageN, posR = pageR*pageN;
//            System.out.println("\t\t\t"+posL+"\t"+posR);

            merge_time-=new Date().getTime();
            int buf_kll_num = 1;
            List<KLLSketchForQuantile> buf_kll_list = new ArrayList<>(buf_kll_num);
            HeapLongStrictKLLSketch merge_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
//            HeapLongKLLSketch merge_worker = new HeapLongKLLSketch(maxMemoryByte);
            for(int i=L;i<Math.min(R,posL);i++)
                merge_worker.update(dataToLong(a[i]));
            for(int i=pageL;i<pageR;i++){
                if(buf_kll_list.size()<buf_kll_num)
                    buf_kll_list.add(KLLArr[i]);
                else{
//                    System.out.println("\t\t\t\t??!!!!!!!!????\t\t");
//                    merge_worker.show();
                    merge_worker.mergeWithTempSpace(buf_kll_list);
//                    merge_worker.show();
                    buf_kll_list.clear();
                    buf_kll_list.add(KLLArr[i]);
//                    System.out.println("\t\t\t\t??!!!!!!!!????\t\t"+Math.abs(merge_worker.getN()/2.0-merge_worker.getApproxRank(doubleToLong(merge_worker.getN()/2.0))));
                }
            }
            if(!buf_kll_list.isEmpty())
                merge_worker.mergeWithTempSpace(buf_kll_list);
            for(int i=Math.max(L,posR);i<R;i++)
                merge_worker.update(dataToLong(a[i]));
            merge_time+=new Date().getTime();


            full_time -= new Date().getTime();
            HeapLongStrictKLLSketch full_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
//            HeapLongKLLSketch full_worker = new HeapLongKLLSketch(maxMemoryByte);
            if(TEST_FULL) {
                for (int i = L; i < R; i++)
                    full_worker.update(dataToLong(a[i]));
            }
            full_time += new Date().getTime();

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R-L);
            Arrays.sort(query_a);

            if(cntTestId==0) {
                merge_worker.show();
                full_worker.show();
            }

            if(merge_worker.cntLevel!=full_worker.cntLevel){
//                System.out.println("\t\t\t\t!!!!!!!!!!!!!!!!!!! merge diff LV....\t"+merge_worker.cntLevel+"\t!=\t"+full_worker.cntLevel+"\t\t\tcntMem="+maxMemoryByte+" "+maxMemoryByte/1024+"KB"+"\ttestID:"+cntTestId+"\t\tL,R:"+L+","+R);
//                merge_worker.show();
//                full_worker.show();
                if(merge_worker.cntLevel>full_worker.cntLevel) {
                    System.out.println("MERGE WORSE");
                    T--;
                    merge_worker.showLevelMaxSize();
                    full_worker.showLevelMaxSize();
                    return;
                }
//                return;
            }
            lv_full+=1.0*full_worker.cntLevel/TEST_CASE;
            lv_merge+=1.0*merge_worker.cntLevel/TEST_CASE;
//            System.out.println("|!!|"+merge_worker.findMinValueWithRank((int)(queryN*0.01)));
            double q_add=0.0001,q_start=q_add,q_end=1-q_add,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
            for(double q=q_start;q<q_end+1e-10;q+=q_add){
                int query_rank = (int)(q*queryN);

                double merge_v = longToResult(merge_worker.findMinValueWithRank(query_rank));
                int merge_delta_rank = getDeltaRank(query_a, queryN, merge_v, query_rank);
                int merge_old_delta = getValueActualRank(query_a,queryN,merge_v)-query_rank;//            merge_delta_rank=merge_old_delta;
                double merge_relative_err = 1.0*merge_delta_rank/(queryN);
                err_merge+=Math.abs(merge_relative_err)/(q_count*TEST_CASE);
                MMP_merge+=1.0 *  Math.abs(Math.abs(merge_old_delta)-Math.abs(merge_delta_rank))  /(queryN)/(q_count*TEST_CASE);

                if(TEST_FULL) {
                    double full_v = longToResult(full_worker.findMinValueWithRank(query_rank));
                    int full_delta_rank = getDeltaRank(query_a, queryN, full_v, query_rank);
                    int full_old_delta = getValueActualRank(query_a,queryN,full_v)-query_rank;//          full_delta_rank=full_old_delta;
                    double full_relative_err = 1.0 * full_delta_rank / (queryN);
                    err_full += Math.abs(full_relative_err) / (q_count * TEST_CASE);
                    MMP_full+=1.0 *  Math.abs(Math.abs(full_old_delta)-Math.abs(full_delta_rank))  /(queryN)/(q_count*TEST_CASE);
//                    System.out.println("?\t\tfull_delta:"+full_delta_rank+"\t\tfull_old_delta:"+full_old_delta);
                }

//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("\t\t"+queryN+"\t"+maxMemoryByte/1024+"\t\t"+err_merge+"\t"+err_full);//        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_mergeBuf-err_full)+"\t"+queryN*(err_merge-err_full));


        System.out.println("\t\t\t\t\t\t\tMMP DEBUG"+"\t"+MMP_merge+"\t"+MMP_full);
        err_result.set(RESULT_LINE,err_result.get(RESULT_LINE).concat("\t\t\t\t\t\t"+maxMemoryByte/1024+"\t\t"+err_merge+(TEST_FULL?("\t"+err_full):"")+"\t\t\t"+lv_merge+(TEST_FULL?("\t"+lv_full):"")));
        if(TEST_FULL)err_result.set(RESULT_LINE,err_result.get(RESULT_LINE).concat("\t\t"+err_full/err_merge));
        time_result.set(RESULT_LINE, time_result.get(RESULT_LINE).concat("\t\t\t\t\t\t" + merge_time+(TEST_FULL?("\t"+full_time):"")));

//        err_result[0][dataType]=err_merge;
//        err_result[1][dataType]=err_full;
    }
    public void show_time_result(){
        System.out.println(time_result);
    }
//    public void show_err_result(){
//        System.out.println(err_result);
//    }
    public static void setTestCase(int tc){TEST_CASE=tc;}


//    public void one_trail(int query_N, int query_Mem,int chunk_sketch_Byte) throws IOException{
//
//        for (int dataType = 1; dataType <= 4; dataType++){
//                testMergeError(0, query_N, query_Mem, chunk_sketch_Byte); // CHECK IT
//            }
////        for(int i=0;i<=4;i++)
////            System.out.print(err_result[0][i]+"\t"+err_result[1][i]+"\t\t\t");
////        System.out.println();
//    }


    public static void main(String[] args) throws IOException{
        long START=new Date().getTime();
        IntervalEvaluatingMergingKLLSketch_222 main;
//        main = new MainForMergeStatErrorKLL();
//        main.prepareA();
//        int tmp_seri = 1<<9,tmp_mem=1<<15;
//        main.prepareKLL(tmp_seri);
////        for(int num=1;num<=8;num++)
////            main.testKLL(num,tmp_mem,tmp_seri);
//        for(int num=1;num<=pageNum;num*=2)
//            main.testKLL(num,tmp_mem,tmp_seri);
//        main.show_time_result();

//        System.out.println("interval query"+"\n");
//        for (int dataType = 1; dataType <= 4; dataType++) {
//            main = new IntervalEvaluatingMergingKLLSketch();
//            main.prepareA(dataType);
//            for (int i : new int[]{/*10000000/*,20000000,30000000,40000000,*/50000000})
//                for (int query_mem : new int[]{/*1024*16,*/1024 * 32, 1024 * 64,/*1024*128/*,1024*256/*,1024*512,1024*1024*/})
//                    for (int chunk_seri : new int[]{/*128,256,512,*/1024/*,2048,4096,8192*/}) {
//                        { // CHECK IT
//                            System.out.println("\n\n\t\tKLL\tdataType:" + dataType + "\t");
//                            System.out.println("\tpageN:" + pageN + "\t|summary|:" + chunk_seri + "\t|Memory|:" + query_mem);
//                            main.testMergeError(0, i, query_mem, chunk_seri);
//                            main.show_time_result();
////                    System.out.println(main.resultText);
//                        }
//                    }
//        }

        System.out.println("interval query"+"\n");
        for (int startType=2,endType=2,dataType = startType; dataType <= endType; dataType++){ // CHECK IT
            main = new IntervalEvaluatingMergingKLLSketch_222();
            main.prepareA(dataType);
            for(int queryN : new int[]{/*10000000,20000000,30000000,40000000,*/40000000/**/})
//                for(int query_mem : new int[]{1024*128,1024*144,1024*152,1024*160,1024*168,1024*184,/*1024*176,*/1024*192,/*1024*208,*//*1024*224,*/1024*200,1024*256,1024*384,1024*480,1024*512,1024*640,1024*1024,1024*2048})
//                for(int query_mem : new int[]{1024*168,1024*184,1024*200})
                for(int query_mem : new int[]{/*1024*144,1024*168,1024*256,1024*480,*/1024*640/*,1024*1024,1024*2048*/})
                    for(int chunk_seri : new int[]{4096}) {
//                        System.out.println("\tpageN:" + pageN + "\t|summary|:" + tmp_seri + "\t|Memory|:" + tmp_mem);
//                        main.prepareWorker(chunk_seri);
                        if(dataType==startType) {
                            err_result.add("N:" + queryN+ ", " + "|M_c|:" + chunk_seri + ", " + "M:\t" + query_mem/1024  + "\t");
                            time_result.add("N:" + queryN+ ", " + "|M_c|:" + chunk_seri + ", " + "M:\t" + query_mem/1024  + "\t");
                        }
                        main.testMergeError(queryN, query_mem, chunk_seri); // CHECK IT
//                        System.out.println("");
                        main.RESULT_LINE++;
//                        main.show_time_result();
                    }
        }
        System.out.println("CHUNK-KLL & KLL\nTEST_CASE="+TEST_CASE);
        System.out.println("\nError rate:");
        for(String s:err_result)
            System.out.println(s);
//        System.out.println("\nQuery Time:");
//        for(String s:time_result)
//            System.out.println(s);
        System.out.println("\t\tALL_TIME:"+(new Date().getTime()-START));
    }
}
