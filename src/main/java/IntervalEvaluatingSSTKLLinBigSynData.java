import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class IntervalEvaluatingSSTKLLinBigSynData {
    int dataType = 1;
    static int LSM_T=30;
    static int pageN = 8192*16;
    static int N = 550000000/pageN*pageN+pageN, pageNum=N/pageN; // CHECK IT
    public static int TEST_CASE=2; // CHECK IT
    boolean TEST_FULL=true; // CHECK IT
    static double[] a,sorted_a;
    static KLLSketchForQuantile[] pageKLL;
    static ObjectArrayList<ObjectArrayList<SST>> lsmNode;
//    String time_result="";
    static ArrayList<String> err_result = new ArrayList<>();
    static ArrayList<String> time_result = new ArrayList<>();
    boolean show_time = false,show_err=true;
    int RESULT_LINE=0;
    Random random = new Random(233);

    public class SST{
        public int level;
        public ObjectArrayList<ObjectArrayList<KLLSketchForQuantile>> sketch;
        public ObjectArrayList<LongArrayList> nodeMinT,nodeMaxT;
        public ObjectArrayList<IntArrayList>nodeStartDataID;
    }

    public void prepareA(int dataType) throws IOException {
        if(a==null){
            a = new double[N];
//            sorted_a = new double[N];
        }
        this.dataType = dataType;

        if(dataType==0){
            for(int i=0;i<N;i++)a[i]=
//                Math.pow(-1,random.nextInt(2))*Math.pow(10.0,(2*Math.pow(random.nextDouble(),2)-1)*300);
                //(double)(i/10000000);
//                random.nextGaussian();
            i;
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
    public void prepareWorker(int maxSeriesByte,int sketchSizeRatio){
//        System.out.println("   prepareWorker start chunk. totalMem:"+1.0*Runtime.getRuntime().totalMemory()/1024/1024 +"\t\tfreeMem:"+ 1.0*Runtime.getRuntime().freeMemory()/1024/1024);
        pageKLL = new KLLSketchForQuantile[pageNum];
        int enoughMemByte = pageN*8;
        for(int i=0;i<pageNum;i++) {
            LongKLLSketch worker = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);
            for (int j = 0; j<pageN; j++) worker.update(dataToLong(a[i*pageN+j]));
            worker.compactBeforeSerialization();
            pageKLL[i] = worker;
//            if(i%100==0)System.out.println("   prepareWorker doing."+1.0*i/pageNum+" totalMem:"+1.0*Runtime.getRuntime().totalMemory()/1024/1024 +"\t\tfreeMem:"+ 1.0*Runtime.getRuntime().freeMemory()/1024/1024);

//            if(i==0)worker.showNum();
        }
//        System.out.println("   prepareWorker start lsm. totalMem:"+1.0*Runtime.getRuntime().totalMemory()/1024/1024 +"\t\tfreeMem:"+ 1.0*Runtime.getRuntime().freeMemory()/1024/1024);

        int sketchNum = 1,SSTlv=0;
        lsmNode = new ObjectArrayList<>();
        lsmNode.add(new ObjectArrayList<>());
        while(sketchNum*LSM_T<=pageNum) {
            sketchNum*=LSM_T;
            SSTlv++;
            lsmNode.add(new ObjectArrayList<>());
        }
//        System.out.println("\t\t??\t"+lsmNode.size());
        for(int remainingPage=pageNum;SSTlv>=0;remainingPage%=sketchNum,SSTlv--,sketchNum/=LSM_T){
            for(int i=0;i<remainingPage/sketchNum;i++){
                SST sst = new SST();
                lsmNode.get(SSTlv).add(sst);
                sst.sketch = new ObjectArrayList<>();
                sst.nodeMinT = new ObjectArrayList<>();
                sst.nodeMaxT = new ObjectArrayList<>();
                sst.nodeStartDataID = new ObjectArrayList<>();

                for(int sketchLV=0,width=sketchNum;sketchLV<=SSTlv;sketchLV++,width/=LSM_T) {
                    sst.sketch.add(new ObjectArrayList<>());
                    sst.nodeMinT.add(new LongArrayList());
                    sst.nodeMaxT.add(new LongArrayList());
                    sst.nodeStartDataID.add(new IntArrayList());
                    if(sketchLV==0) {
                        for (int j = 0; j < width; j++) {
                            int pageID = pageNum - remainingPage + i * sketchNum + j;
                            sst.sketch.get(0).add(pageKLL[pageID]);
                            sst.nodeMinT.get(0).add(pageID * pageN);
                            sst.nodeMaxT.get(0).add(pageID * pageN + pageN - 1);
                            sst.nodeStartDataID.get(0).add(pageID * pageN);
                        }
                    }else{
                        for(int j=0;j<width;j++){
                            KLLSketchForSST sketch = new KLLSketchForSST();
                            for(int k=0;k<LSM_T;k++) {
                                sketch.addSubSketch(sst.sketch.get(sketchLV - 1).get(j * LSM_T + k));
                            }
                            sketch.compactSubSketches(sketchSizeRatio);
                            sst.sketch.get(sketchLV).add(sketch);
                            long mn = sst.nodeMinT.get(sketchLV-1).getLong(j * LSM_T);
                            long mx = sst.nodeMaxT.get(sketchLV-1).getLong(j * LSM_T+LSM_T-1);
                            sst.nodeMinT.get(sketchLV).add(mn);
                            sst.nodeMaxT.get(sketchLV).add(mx);
                            sst.nodeStartDataID.get(sketchLV).add(sst.nodeStartDataID.get(sketchLV-1).getInt(j * LSM_T));
//                            System.out.println("\t\ta sketch in sst.  T:"+mn+"..."+mx+"\t\t");
//                            System.out.print("\t\t\t\t");sketch.show();
                        }
                    }
                }
            }
        }
//        for(int lv=lsmNode.size()-1;lv>=0;lv--) {
//            lsmNode.get(lsmNode.size() - 1).get(0).sketch.get(lv).get(0).show();
//        }
//        for(int lsmLV=lsmNode.size()-1;lsmLV>=0;lsmLV--) {
//            System.out.println("\t\t\t\tthere are "+lsmNode.get(lsmLV).size()+" lv "+lsmLV+" lsm components.");
//        }
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

    private boolean inInterval(long x, long y, long L, long R) {
        return x >= L && y <= R;
    }

    private boolean inInterval(long x, long L, long R) {
        return x >= L && x <= R;
    }

    private boolean overlapInterval(long x, long y, long L, long R) { // [L,R]
        return !(y < L || x > R);
    }

    private void range_query(int lv, int p, SST sst, HeapLongStrictKLLSketch query_sketch, long L, long R, long otherL, long otherR){ // []

        long cntL=sst.nodeMinT.get(lv).getLong(p),cntR=sst.nodeMaxT.get(lv).getLong(p);
        if(!overlapInterval(cntL,cntR,L,R))return;
//            System.out.println("\t\t\t\t\tcntSST T:"+cntL+"..."+cntR+"\t\t\t\tqueryLR:"+L+","+R);
        if(inInterval(cntL,cntR,L,R)&&inInterval(cntL,cntR,otherL,otherR)){
//            System.out.println("\t\tmerge with T:"+cntL+"..."+cntR+"\t\t\tlv="+lv+"\t\tcntN:"+query_sketch.getN());
            query_sketch.mergeWithTempSpace(Collections.singletonList(sst.sketch.get(lv).get(p)));
//            query_sketch.show();
            return;
        }
        if(lv>0) {
            for(int i=0;i<LSM_T;i++)
                range_query(lv - 1, p*LSM_T+i, sst, query_sketch, L, R, otherL, otherR);
        }else{
            for (int start=sst.nodeStartDataID.get(lv).getInt(p),j = start; j < start + pageN; j++) {
                int index = j;//cc.getInt(j);
                if (inInterval(index, L, R))
                    query_sketch.update(dataToLong(a[index]));
            }
        }
    }



    private long dataToLong(double data)  {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }

    private double longToResult(long result)  {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }

    public void testMergeError(int queryN, int maxMemoryByte, int maxSeriesByte, int sketchSizeRatio) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time=0, merge_sst_time=0;
        double err_full=0,err_merge_sst=0;
        double lv_full=0,lv_merge=0;
        double MMP_merge = 0,MMP_full = 0;

        sorted_a = null;
        System.gc();
        try {
            Thread.sleep(2000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        sorted_a = new double[queryN];

        int[] LL = new int[TEST_CASE*10];
        int[] RR = new int[TEST_CASE*10];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE*10; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }

        for(int T=0,cntTestId=0;T<TEST_CASE;T++,cntTestId++){
//            if((cntTestId&1)==0)//CHECK IT!
                prepareWorker(maxSeriesByte,sketchSizeRatio);
            int L = LL[cntTestId], R = RR[cntTestId];
            int pageL = (L+pageN-1)/pageN, pageR = R/pageN;
            int posL = pageL*pageN, posR = pageR*pageN;
//            System.out.println("\t\t\t"+posL+"\t"+posR);


            merge_sst_time-=new Date().getTime();
            HeapLongStrictKLLSketch merge_sst_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            for(int lsmLV=lsmNode.size()-1;lsmLV>=0;lsmLV--){
                for(SST sst:lsmNode.get(lsmLV))
                    range_query(lsmLV,0,sst,merge_sst_worker,L,R-1,Long.MIN_VALUE,Long.MAX_VALUE);
//                System.out.println("--------------------------------------------lsm LV over. cntN:"+merge_sst_worker.getN());
            }
            merge_sst_time+=new Date().getTime();

            full_time -= new Date().getTime();
            HeapLongStrictKLLSketch full_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            for (int i = L; i < R; i++)
                full_worker.update(dataToLong(a[i]));
            full_time += new Date().getTime();

//            System.out.println("????????????????????? "+merge_sst_worker.getN()+"\t\t"+full_worker.getN());
            assert merge_sst_worker.getN()==full_worker.getN();


            if (R - L >= 0) System.arraycopy(a, L, sorted_a, 0, R-L);
            Arrays.sort(sorted_a);

//            if(cntTestId==0) {
//                merge_sst_worker.show();
//                full_worker.show();
//            }
//
//            lv_full+=1.0*full_worker.cntLevel/TEST_CASE;
            lv_merge+=1.0*merge_sst_worker.cntLevel/TEST_CASE;
//            System.out.println("|!!|"+merge_worker.findMinValueWithRank((int)(queryN*0.01)));
            double q_add=0.0001,q_start=q_add,q_end=1-q_add,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
            for(double q=q_start;q<q_end+1e-10;q+=q_add){
                int query_rank = (int)(q*queryN);

                double merge_v = longToResult(merge_sst_worker.findMinValueWithRank(query_rank));
                int merge_delta_rank = getDeltaRank(sorted_a, queryN, merge_v, query_rank);
                double merge_relative_err = 1.0*merge_delta_rank/(queryN);
                err_merge_sst+=Math.abs(merge_relative_err)/(q_count*TEST_CASE);

                if(TEST_FULL) {
                    double full_v = longToResult(full_worker.findMinValueWithRank(query_rank));
                    int full_delta_rank = getDeltaRank(sorted_a, queryN, full_v, query_rank);
                    double full_relative_err = 1.0 * full_delta_rank / (queryN);
                    err_full += Math.abs(full_relative_err) / (q_count * TEST_CASE);
//                    System.out.println("?\t\tfull_delta:"+full_delta_rank+"\t\tfull_old_delta:"+full_old_delta);
                }

//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
//        System.out.println("\t\t"+queryN+"\t"+maxMemoryByte/1024+"\tT_s\t"+sketchSizeRatio+"\t\t"+err_merge_sst+"\t"+err_full);//        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_mergeBuf-err_full)+"\t"+queryN*(err_merge-err_full));


        err_result.set(RESULT_LINE,err_result.get(RESULT_LINE).concat("\t\t\t\t\t\t"+err_merge_sst+(TEST_FULL?("\t"+err_full):"")+"\t\t\t"+lv_merge+(TEST_FULL?("\t"+lv_full):"")));
        if(TEST_FULL)err_result.set(RESULT_LINE,err_result.get(RESULT_LINE).concat("\t\t"+err_full/err_merge_sst));
        time_result.set(RESULT_LINE, time_result.get(RESULT_LINE).concat("\t\t\t\t\t\t" + merge_sst_time+(TEST_FULL?("\t"+full_time):"")));

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
        IntervalEvaluatingSSTKLLinBigSynData main;
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
        for (int startType=0,endType=0,dataType = startType; dataType <= endType; dataType++){ // CHECK IT
            main = new IntervalEvaluatingSSTKLLinBigSynData();
            main.prepareA(dataType);
            for(int queryN : new int[]{100000000,200000000,300000000,400000000,500000000})
//                for(int query_mem : new int[]{1024*128,1024*144,1024*152,1024*160,1024*168,1024*184,/*1024*176,*/1024*192,/*1024*208,*//*1024*224,*/1024*200,1024*256,1024*384,1024*480,1024*512,1024*640,1024*1024,1024*2048})
//                for(int query_mem : new int[]{1024*168,1024*184,1024*200})
                for(int query_mem : new int[]{/*1024*144,1024*168,*/8192*32*16/*,1024*480,1024*640,1024*1024,1024*2048/*1024*512*/})
                    for(int chunk_seri : new int[]{4096*16})
                        for(int SketchSizeRatio:new int[]{/*1,2,*/4/*,8,16,32/*8/**/})
                        {
//                        System.out.println("\tpageN:" + pageN + "\t|summary|:" + tmp_seri + "\t|Memory|:" + tmp_mem);
//                        main.prepareWorker(chunk_seri);
                        if(dataType==startType) {
                            err_result.add("N:" + queryN+ ", " + "T_s:" + SketchSizeRatio + ", " + "M:\t" + query_mem/1024  + "\t");
                            time_result.add("N:" + queryN+ ", " + "T_s:" + SketchSizeRatio + ", " + "M:\t" + query_mem/1024  + "\t");
                        }
//                        main.TEST_FULL=SketchSizeRatio==1;
                        main.testMergeError(queryN, query_mem, chunk_seri, SketchSizeRatio); // CHECK IT
//                        System.out.println("");
                        main.RESULT_LINE++;
//                        main.show_time_result();
                    }
        }
        System.out.println("SST-KLL & Online-KLL\nTEST_CASE="+TEST_CASE);
        System.out.println("\nError rate:\t\t\t\t\t\t\t\t\t\t\t\tsst-kll err\t\t\tonline-kll err\t\t\t\tsst-avg_lv\t\t\tonline-avg_lv\t\t\tonlineErr/sstErr");
        for(String s:err_result)
            System.out.println(s);
//        System.out.println("\nQuery Time:");
//        for(String s:time_result)
//            System.out.println(s);
        System.out.println("\t\tALL_TIME:"+(new Date().getTime()-START));
    }
}
