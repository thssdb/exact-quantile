import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;

public class ComparingBuildingOptimalAndOriginalSketch {
    int dataType;
    static int pageN = 8192;
    static int N = 55000000/pageN*pageN, pageNum=N/pageN; // CHECK IT
    public static int TEST_CASE=1; // CHECK IT
    static double[] a;
    static KLLSketchForQuantile[] OptimalArr;
    static HeapLongStrictKLLSketch[] OriginalArr;
    static ArrayList<String> err_result = new ArrayList<>();
    static ArrayList<String> time_result = new ArrayList<>();
    static ArrayList<String> build_time_result = new ArrayList<>();
    boolean show_time = false,show_err=false;
    int RESULT_LINE=0;

    Random random = new Random(233);

    public void prepareA(int dataType) throws IOException {
        if(a==null)a = new double[N];
        this.dataType = dataType;

        if(dataType==0){
            for(int i=0;i<N;i++)a[i]=Math.pow(-1,random.nextInt(2))*Math.pow(10.0,(2*Math.pow(random.nextDouble(),2)-1)*300);
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
    public void prepareWorker(int maxSeriesByte,int TEST_CASE){
        long original_time=0, optimal_time=0;
        int enoughMemByte = pageN*9;

        int byteForOriginal = maxSeriesByte;
//        for(int cntByte = maxSeriesByte;cntByte<=maxSeriesByte*8;cntByte+=8) {
//            HeapLongStrictKLLSketch original_worker = new HeapLongStrictKLLSketch(cntByte);
//            for (int j = 0; j < pageN; j++) original_worker.update(dataToLong(j));
//            if(original_worker.getNumLen()*8<=maxSeriesByte)
//                byteForOriginal = cntByte;
//        }
//        System.out.println("\t\t\t\t\tbyteForOriginal:"+byteForOriginal);

        for(int T=0;T<TEST_CASE;T++) {
            OptimalArr = new KLLSketchForQuantile[pageNum];
            OriginalArr = new HeapLongStrictKLLSketch[pageNum];
            for (int i = 0; i < pageNum; i++) {
                OptimalArr[i] = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);
                OriginalArr[i] = new HeapLongStrictKLLSketch(byteForOriginal);
            }

            optimal_time -= new Date().getTime();
            for (int i = 0; i < pageNum; i++) {
//            LongKLLSketch worker = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);
                for (int j = 0; j < pageN; j++) OptimalArr[i].update(dataToLong(a[i * pageN + j]));
                ((LongKLLSketch) OptimalArr[i]).compactBeforeSerialization();
//            OptimalArr[i] = worker;
//            if(i==0){
//                worker.show();
//                worker = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);
//                for (int j = 0; j<pageN; j++) worker.update(j);
//                worker.compactBeforeSerialization();
//                worker.showNum();
//            }
            }
            optimal_time += new Date().getTime();

            original_time -= new Date().getTime();
            for (int i = 0; i < pageNum; i++) {
//            if(i==0) {
//            }
//            HeapLongStrictKLLSketch original_worker = new HeapLongStrictKLLSketch(byteForOriginal);
                for (int j = 0; j < pageN; j++) OriginalArr[i].update(dataToLong(a[i * pageN + j]));
//            OriginalArr[i] = original_worker;
//            if(i==0){
//                original_worker.show();
//                original_worker = new HeapLongStrictKLLSketch(byteForOriginal);
//                for (int j = 0; j < pageN; j++) original_worker.update(j);
//                original_worker.showNum();
//            }
            }
            original_time += new Date().getTime();
        }
        if (show_time) {
//            System.out.print("\t\t"+(new Date().getTime()-ST));
            build_time_result.set(RESULT_LINE, build_time_result.get(RESULT_LINE).concat("\t\t\t" + (optimal_time/TEST_CASE)));
            build_time_result.set(RESULT_LINE, build_time_result.get(RESULT_LINE).concat("\t" + (original_time/TEST_CASE)));
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
        double original_time=0, optimal_time=0;
        double err_original=0,err_optimal=0;
        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = 0;//random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }

//        OptimalArr[0].show();
//        OriginalArr[0].show();

        for(int T=0;T<TEST_CASE;T++){
            if((T&1)==0) {
                show_time=false;
                prepareWorker(maxSeriesByte, 1);
            }
            int L = LL[T], R = RR[T];
            int pageL = (L+pageN-1)/pageN, pageR = R/pageN;
            int posL = pageL*pageN, posR = pageR*pageN;
            int merge_buf_size = maxMemoryByte/32/maxSeriesByte;
            ObjectArrayList<KLLSketchForQuantile> merge_buf=new ObjectArrayList<>(merge_buf_size);
//            System.out.println("\t\t\t"+posL+"\t"+posR);

            HeapLongStrictKLLSketch optimal_worker=null;
//            optimal_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
//            optimal_time-=new Date().getTime();
//            for(int i=L;i<Math.min(R,posL);i++)
//                optimal_worker.update(dataToLong(a[i]));
//            for(int i=pageL;i<pageR;/*i++*/) {
////                optimal_worker.mergeWithTempSpace(OptimalArr[i]);
//                merge_buf.clear();
//                for(int j=0;j<merge_buf_size&&i+j<pageR;j++)merge_buf.add(OptimalArr[i+j]);
//                optimal_worker.mergeWithTempSpace(merge_buf);
//                i+=merge_buf_size;
//            }
//            for(int i=Math.max(L,posR);i<R;i++)
//                optimal_worker.update(dataToLong(a[i]));
//            optimal_time+=new Date().getTime();

            for(int merge_t=0;merge_t<TEST_CASE;merge_t++){
                optimal_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
                optimal_time-=1.0*new Date().getTime()/TEST_CASE;
                for(int i=L;i<Math.min(R,posL);i++)
                    optimal_worker.update(dataToLong(a[i]));
//                for(int i=pageL;i<pageR;i++)
//                    optimal_worker.mergeWithTempSpace(OptimalArr[i]);
                for(int i=pageL;i<pageR;i+=merge_buf_size/*i++*/) {
//                optimal_worker.mergeWithTempSpace(OptimalArr[i]);
                    merge_buf.clear();
                    for(int j=0;j<merge_buf_size&&i+j<pageR;j++)merge_buf.add(OptimalArr[i+j]);
                    optimal_worker.mergeWithTempSpace(merge_buf);
                }
                for(int i=Math.max(L,posR);i<R;i++)
                    optimal_worker.update(dataToLong(a[i]));
                optimal_time+=1.0*new Date().getTime()/TEST_CASE;
            }



            HeapLongStrictKLLSketch original_worker=null;
//            original_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
//            original_time-=new Date().getTime();
//            for(int i=L;i<Math.min(R,posL);i++)
//                original_worker.update(dataToLong(a[i]));
//            for(int i=pageL;i<pageR;i+=merge_buf_size/*i++*/) {
////                optimal_worker.mergeWithTempSpace(OptimalArr[i]);
//                merge_buf.clear();
//                for(int j=0;j<merge_buf_size&&i+j<pageR;j++)merge_buf.add(OriginalArr[i+j]);
//                original_worker.mergeWithTempSpace(merge_buf);
//            }
//            for(int i=Math.max(L,posR);i<R;i++)
//                original_worker.update(dataToLong(a[i]));
//            original_time+=new Date().getTime();

            for(int merge_t=0;merge_t<TEST_CASE;merge_t++){
                original_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
                original_time-=1.0*new Date().getTime()/TEST_CASE;
                for(int i=L;i<Math.min(R,posL);i++)
                    original_worker.update(dataToLong(a[i]));
                for(int i=pageL;i<pageR;i+=merge_buf_size/*i++*/) {
//                optimal_worker.mergeWithTempSpace(OptimalArr[i]);
                    merge_buf.clear();
                    for(int j=0;j<merge_buf_size&&i+j<pageR;j++)merge_buf.add(OriginalArr[i+j]);
                    original_worker.mergeWithTempSpace(merge_buf);
                }
                for(int i=Math.max(L,posR);i<R;i++)
                    original_worker.update(dataToLong(a[i]));
                original_time+=1.0*new Date().getTime()/TEST_CASE;
            }

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R-L);
            Arrays.sort(query_a);

            double q_add=0.0001,q_start=q_add,q_end=1-q_add,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
            for(double q=q_start;q<q_end+1e-10;q+=q_add){
                int query_rank = (int)(q*queryN);

                double optimal_v = longToResult(optimal_worker.findMinValueWithRank(query_rank));
                int optimal_delta_rank = getDeltaRank(query_a, queryN, optimal_v, query_rank);
                double optimal_relative_err = 1.0*optimal_delta_rank/(queryN);
                err_optimal+=Math.abs(optimal_relative_err)/(q_count*TEST_CASE);

                double original_v = longToResult(original_worker.findMinValueWithRank(query_rank));
                int original_delta_rank = getDeltaRank(query_a, queryN, original_v, query_rank);
                double original_relative_err = 1.0*original_delta_rank/(queryN);
                err_original+=Math.abs(original_relative_err)/(q_count*TEST_CASE);

//                System.out.println("?\t\toriginal:"+original_v+" delta:"+original_delta_rank+"\t\toptimal:"+optimal_v+" delta:"+optimal_delta_rank);
            }
        }
        if(show_err) {
            System.out.println("\t\t\t" + err_optimal + "\t" + err_original);
        }
        err_result.set(RESULT_LINE,err_result.get(RESULT_LINE).concat("\t\t\t"+err_optimal+"\t"+err_original));
        if(TEST_CASE!=0&&show_time) {
//            if(show_time) {
//                System.out.print("\t" + optimal_time / TEST_CASE);
//                System.out.print("\t" + original_time / TEST_CASE);
//                System.out.println();
//            }
            time_result.set(RESULT_LINE, time_result.get(RESULT_LINE).concat("\t\t\t" + optimal_time / TEST_CASE+"\t"+original_time / TEST_CASE));
        }
        //        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_optimalBuf-err_original)+"\t"+queryN*(err_optimal-err_original));

    }
    public void show_time_result(){
        System.out.println(time_result);
    }
    public static void setTestCase(int tc){TEST_CASE=tc;}




    public static void main(String[] args) throws IOException{
        long START=new Date().getTime();
        ComparingBuildingOptimalAndOriginalSketch main;
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
        for (int startType=1,endType=1,dataType = startType; dataType <= endType; dataType++){ // CHECK IT
            main = new ComparingBuildingOptimalAndOriginalSketch();
            main.prepareA(dataType);
            for(int queryN : new int[]{/*10000000,20000000,30000000,40000000,*/40000000})
                for(int query_mem : new int[]{/*1024*16,1024*32,1024*64,*/1024*256/*,1024*256*//*,1024*512,1024*1024*/}) {
//                    main.show_time = false;
//                    for (int chunk_seri : new int[]{512,1024,2048,4096,8192})
//                        main.prepareWorker(chunk_seri); // for time test.
//                    main.show_time = true;
                    for (int chunk_seri : new int[]{256/*,512,1024,2048,4096*/}) {
                        if (dataType == startType) {
                            err_result.add("N:" + queryN + ", " + "M:" + query_mem + ", " + "|M_c|:" + chunk_seri + "\t");
                            time_result.add("N:" + queryN + ", " + "M:" + query_mem + ", " + "|M_c|:" + chunk_seri + "\t");
                            build_time_result.add("N:" + queryN + ", " + "M:" + query_mem + ", " + "|M_c|:" + chunk_seri + "\t");
                        }
                        main.show_time = false;main.prepareWorker(chunk_seri,1);
                        main.show_time = true;main.prepareWorker(chunk_seri,TEST_CASE);
                        main.testMergeError(queryN, query_mem, chunk_seri); // CHECK IT
//                        System.out.println("");
                        main.RESULT_LINE++;
//                        main.show_time_result();
                    }
                }
            }
        System.out.println("Comparing Chunk Sketch: optimal & original\nTEST_CASE="+TEST_CASE);
        System.out.println("\nTime cost for construction:");
        for(String s:build_time_result)
            System.out.println(s);
        System.out.println("\nError Rate of Merge Result:");
        for(String s:err_result)
            System.out.println(s);
        System.out.println("\nTime cost For Merging pre-computed sketches:");
        for(String s:time_result)
            System.out.println(s);
        System.out.println("\t\tALL_TIME:"+(new Date().getTime()-START));
    }
}
