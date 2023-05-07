import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntDoubleImmutablePair;
import it.unimi.dsi.fastutil.ints.IntDoublePair;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class IntervalEvaluatingLatencySSTKLL {
    int dataType = 1;
    static int LSM_T=30;
    static int pageN = 8192;
//    static double[] muList = new double[]{2.0};
//    static double[] sigList = new double[]{0.0,2.0,2.5,2.6,2.7,2.9,3.0,3.2};
    static double[] muList = new double[]{3.0};
//    static double[] sigList = new double[]{0.0,1.0,1.5,2.0,2.2,2.4,2.5,2.6};
static double[] sigList = new double[]{/*0.0,1.0,2.1,2.4,*/2.6};
//    static double[] sigList = new double[]{/*0.0,1.0,2.1,2.4,*/2.6};
    static int N = 55000000/pageN*pageN, pageNum=N/pageN; // CHECK IT
    public static int TEST_CASE=64; // CHECK IT
    boolean TEST_FULL = false; // CHECK IT
    static double[] a;
    static KLLSketchForQuantile[] pageKLL;
    static long[] workerMinT, workerMaxT;
    static ObjectArrayList<ObjectArrayList<SST>> lsmNode;
    //    String time_result="";
    static ArrayList<String> err_result = new ArrayList<>();
    static ArrayList<String> time_result = new ArrayList<>();
    boolean show_time = false, show_err = true;
    int RESULT_LINE = 0;
    Random random = new Random(2333);
    static long[] aa;
    static long[] bb;
    static IntArrayList cc;
    int mergedSketchN;


    public class SST{
        public int level;
        public ObjectArrayList<ObjectArrayList<KLLSketchForQuantile>> sketch;
        public ObjectArrayList<LongArrayList> nodeMinT,nodeMaxT;
        public ObjectArrayList<IntDoublePair>compactedData;
    }

    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] =
                i;
//                    Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
            return;
        }
        BufferedReader reader = null;
        if (dataType == 1)
            reader = new BufferedReader(new FileReader(new File("1_bitcoin.csv")));
        if (dataType == 2)
            reader = new BufferedReader(new FileReader(new File("2_SpacecraftThruster.txt")));
        if (dataType == 3)
            reader = new BufferedReader(new FileReader(new File("3_taxipredition8M.txt")));
        if (dataType == 4)
            reader = new BufferedReader(new FileReader(new File("4_wh.csv")));
        assert reader != null;
        reader.readLine(); // ignore first line.
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if (cntN == N) break;
        }
    }

    public void prepareDisorder(double mu,double sig) {
        random = new Random(233);
        aa = new long[N];
        bb = new long[N];
        cc = new IntArrayList(N);
        for (int i = 0; i < N; i++) {
            cc.add(i);
            aa[i] = i;
            bb[i] = (long) Math.round(i + Math.exp(mu + sig * random.nextGaussian()));
        }
        cc.sort((x, y) -> (Long.compare(bb[x], bb[y])));
    }

    public void prepareWorker(int maxSeriesByte,int sketchSizeRatio) {
        int enoughMemByte = pageN * 9;
//        pageKLL = new KLLSketchForQuantile[pageNum];
//        workerMinT = new long[pageNum];
//        workerMaxT = new long[pageNum];
//        for (int i = 0; i < pageNum; i++) {
//            LongKLLSketch worker = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);
//
//            workerMinT[i] = Long.MAX_VALUE;
//            workerMaxT[i] = Long.MIN_VALUE;
//            for (int j = 0; j < pageN; j++) {
//                int index = cc.getInt(i * pageN + j);
//                worker.update(dataToLong(a[index]));
//                workerMinT[i] = Math.min(workerMinT[i], index);
//                workerMaxT[i] = Math.max(workerMaxT[i], index);
//            }
//            worker.compactBeforeSerialization();
//            pageKLL[i] = worker;
////            if(i==0)worker.show();
//        }

        int sketchNum = 1,SSTlv=0;
        lsmNode = new ObjectArrayList<>();
        lsmNode.add(new ObjectArrayList<>());
        while(sketchNum*LSM_T<=pageNum) {
            sketchNum*=LSM_T;
            SSTlv++;
            lsmNode.add(new ObjectArrayList<>());
        }
//        System.out.println("\t\t??\t"+lsmNode.size());
        int cntCompactedPos=0;
        for(int lsmLV=SSTlv;lsmLV>=0;lsmLV--){
            int chunkNum = (int)Math.pow(LSM_T,lsmLV);
            while(cntCompactedPos+chunkNum<=pageNum){
                SST sst = new SST();
                lsmNode.get(lsmLV).add(sst);
                ObjectArrayList<IntDoublePair>compactedData =new ObjectArrayList<>();
                for(int i=0;i<chunkNum*pageN;i++){
                    int pos =cntCompactedPos*pageN+i;
                    IntDoublePair pair = new IntDoubleImmutablePair(cc.getInt(pos),a[cc.getInt(pos)]);
                    compactedData.add(pair);
//                    if(compactedData.size()%100==0)System.out.println("\t\t\t\tdata:"+pair.firstInt()+"\t\tv="+pair.secondDouble());
                }
//                System.out.println("\t\t\t\t\tre-sorted data   range:\t"+cntCompactedPos*pageN+"......"+(cntCompactedPos*pageN+chunkNum*pageN));
                compactedData.sort(Comparator.comparingInt(IntDoublePair::firstInt));
                sst.compactedData  = compactedData;
//                System.out.println("\t\t\t\t\t\tsst mnTmxT:"+compactedData.get(0).firstInt()+"..."+compactedData.get(compactedData.size()-1).firstInt());

                sst.sketch = new ObjectArrayList<>();
                sst.nodeMinT = new ObjectArrayList<>();
                sst.nodeMaxT = new ObjectArrayList<>();

                sst.sketch.add(new ObjectArrayList<>());
                sst.nodeMinT.add(new LongArrayList());
                sst.nodeMaxT.add(new LongArrayList());
                for(int i=0;i<chunkNum;i++){
                    LongKLLSketch worker = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);
                    long minT=compactedData.get(i*pageN).firstInt();
                    long maxT=compactedData.get(i*pageN+pageN-1).firstInt();
                    for(int j=0;j<pageN;j++){
                        int t=compactedData.get(i*pageN+j).firstInt();
                        double v = compactedData.get(i*pageN+j).secondDouble();
                        worker.update(dataToLong(v));
                    }
                    worker.compactBeforeSerialization();
//                    int pageID = cntCompactedPos+i;
                    sst.sketch.get(0).add(worker);
                    sst.nodeMinT.get(0).add(minT);
                    sst.nodeMaxT.get(0).add(maxT);
                }
                //System.out.println("\t\t\t\tcnt SST   bottom chunkNum="+chunkNum);
                for(int lv=1,width=chunkNum/LSM_T;lv<=lsmLV;lv++,width/=LSM_T){
                    sst.sketch.add(new ObjectArrayList<>());
                    sst.nodeMinT.add(new LongArrayList());
                    sst.nodeMaxT.add(new LongArrayList());
                    for(int j=0;j<width;j++){
                        KLLSketchForSST sketch = new KLLSketchForSST();
                        for(int k=0;k<LSM_T;k++) {
                            sketch.addSubSketch(sst.sketch.get(lv - 1).get(j * LSM_T + k));
                        }
//                        sketch.show();
//                        sketch.showNum();
                        sketch.compactSubSketches(sketchSizeRatio);
                        sst.sketch.get(lv).add(sketch);
                        long mn = sst.nodeMinT.get(lv-1).getLong(j * LSM_T);
                        long mx = sst.nodeMaxT.get(lv-1).getLong(j * LSM_T+LSM_T-1);
                        sst.nodeMinT.get(lv).add(mn);
                        sst.nodeMaxT.get(lv).add(mx);
//                        sst.nodeStartDataID.get(lv).add(sst.nodeStartDataID.get(lv-1).getInt(j * LSM_T));

//                        if(lv==lsmLV) {
//                            System.out.println("\t\ta sketch in sst.  T:" + mn + "..." + mx + "\t\t");
//                            System.out.print("\t\t\t\t");
//                            sketch.show();
//                            sketch.showNum();
//                        }
                    }
                }

                cntCompactedPos+=chunkNum;
            }
        }
//        for(int remainingPage=pageNum;SSTlv>=0;remainingPage%=sketchNum,SSTlv--,sketchNum/=LSM_T){
//            for(int i=0;i<remainingPage/sketchNum;i++){
//                SST sst = new SST();
//                lsmNode.get(SSTlv).add(sst);
//                sst.sketch = new ObjectArrayList<>();
//                sst.nodeMinT = new ObjectArrayList<>();
//                sst.nodeMaxT = new ObjectArrayList<>();
//                sst.nodeStartDataID = new ObjectArrayList<>();
//
//                for(int sketchLV=0,width=sketchNum;sketchLV<=SSTlv;sketchLV++,width/=LSM_T) {
//                    sst.sketch.add(new ObjectArrayList<>());
//                    sst.nodeMinT.add(new LongArrayList());
//                    sst.nodeMaxT.add(new LongArrayList());
//                    sst.nodeStartDataID.add(new IntArrayList());
//                    if(sketchLV==0) {
//                        for (int j = 0; j < width; j++) {
//                            int pageID = pageNum - remainingPage + i * sketchNum + j;
//                            sst.sketch.get(0).add(pageKLL[pageID]);
//                            sst.nodeMinT.get(0).add(workerMinT[pageID]);
//                            sst.nodeMaxT.get(0).add(workerMaxT[pageID]);
//                            sst.nodeStartDataID.get(0).add(pageID * pageN);
//                        }
//                    }else{
//                        for(int j=0;j<width;j++){
//                            KLLSketchForSST sketch = new KLLSketchForSST();
//                            for(int k=0;k<LSM_T;k++) {
//                                sketch.addSubSketch(sst.sketch.get(sketchLV - 1).get(j * LSM_T + k));
//                            }
//                            sketch.compactSubSketches(sketchSizeRatio);
//                            sst.sketch.get(sketchLV).add(sketch);
//                            long mn = sst.nodeMinT.get(sketchLV-1).getLong(j * LSM_T);
//                            long mx = sst.nodeMaxT.get(sketchLV-1).getLong(j * LSM_T+LSM_T-1);
//                            sst.nodeMinT.get(sketchLV).add(mn);
//                            sst.nodeMaxT.get(sketchLV).add(mx);
//                            sst.nodeStartDataID.get(sketchLV).add(sst.nodeStartDataID.get(sketchLV-1).getInt(j * LSM_T));
//                            System.out.println("\t\ta sketch in sst.  T:"+mn+"..."+mx+"\t\t");
//                            System.out.print("\t\t\t\t");sketch.show();
//                        }
//                    }
//                }
//            }
//        }

//        for(int i=0;i<100;i++)System.out.println(workerMinT[i]+"---"+workerMaxT[i]);
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


    private long dataToLong(double data) {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }

    private double longToResult(long result) {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
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

    private void range_query(int lv, int p, SST sst, HeapLongStrictKLLSketch query_sketch, long L, long R, ObjectArrayList<SST>otherFiles){ // []

        long cntL=sst.nodeMinT.get(lv).getLong(p),cntR=sst.nodeMaxT.get(lv).getLong(p);
//        System.out.println("\t\t\t\t\tcntSST T:"+cntL+"..."+cntR+"\t\t\t\tqueryLR:"+L+","+R);
        if(!inInterval(cntL,cntR,L,R)&&!overlapInterval(cntL,cntR,L,R))return;
//            System.out.println("\t\t\t\t\t!!cntSST T:"+cntL+"..."+cntR+"\t\t\t\tqueryLR:"+L+","+R);
        if(inInterval(cntL,cntR,L,R)){
            boolean overlapped = false;
            for(SST otherFile:otherFiles)
                overlapped|=overlapInterval(cntL,cntR,otherFile.nodeMinT.get(otherFile.sketch.size()-1).getLong(0),otherFile.nodeMaxT.get(otherFile.sketch.size()-1).getLong(0));
            if(!overlapped) {
                mergedSketchN+=sst.sketch.get(lv).get(p).getN();
//                System.out.println("\t\tmerge with T:"+cntL+"..."+cntR+"\t\t\tlv="+lv+"\t\tcntQueryN:"+query_sketch.getN());
                query_sketch.mergeWithTempSpace(Collections.singletonList(sst.sketch.get(lv).get(p)));
//                query_sketch.show();
                return;
            }
        }
        if(lv>0) {
            for(int i=0;i<LSM_T;i++)
                range_query(lv - 1, p*LSM_T+i, sst, query_sketch, L, R, otherFiles);
        }else{
            for(int i=0;i<pageN;i++){
                int t=sst.compactedData.get(p*pageN+i).firstInt();
                if(inInterval(t,L,R)) {
                    query_sketch.update(dataToLong(sst.compactedData.get(p * pageN + i).secondDouble()));
//                    if(query_sketch.getN()%100==0)query_sketch.showNum();
                }
            }
//            for (int start=sst.nodeStartDataID.get(lv).getInt(p),j = start; j < start + pageN; j++) {
//                int index = j;//cc.getInt(j);
//                if (inInterval(index, L, R))
//                    query_sketch.update(dataToLong(a[index]));
//            }
        }
    }


    public void testMergeError(int queryN, int maxMemoryByte, int maxSeriesByte,double sig,int sketchSizeRate) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time = 0, merge_time = 0;
        double err_merge_sst = 0, err_merge = 0;
        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        double  avgAvailable=0;
        mergedSketchN=0;
        for (int T = 0; T < TEST_CASE; T++) {
            if((T&1)==0){
                prepareWorker(maxSeriesByte,sketchSizeRate);
            }
            int L = LL[T], R = RR[T];

            HeapLongStrictKLLSketch merge_sst_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            ObjectArrayList<SST> allFiles =new ObjectArrayList<>();
            for(int lsmLV=lsmNode.size()-1;lsmLV>=0;lsmLV--)
                for(SST sst:lsmNode.get(lsmLV))
                    allFiles.add(sst);
            for(int lsmLV=lsmNode.size()-1;lsmLV>=0;lsmLV--){
                for(SST sst:lsmNode.get(lsmLV)) {
                    allFiles.remove(sst);
//                    System.out.println("\t\t\t\t\tlsm divide.  lsmLV="+lsmLV);
                    range_query(lsmLV, 0, sst, merge_sst_worker, L, R - 1, allFiles);
                    allFiles.add(sst);
                }
//                System.out.println("--------------------------------------------lsm LV over. cntN:"+merge_sst_worker.getN());
            }
//            System.out.println("\t\t\tquery l,r:"+L+"\t\t"+R);
//            merge_sst_worker.show();
//            merge_sst_worker.showNum();

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);

            double q_add=0.0001,q_start=q_add,q_end=1-q_add,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
            for(double q=q_start;q<q_end+1e-10;q+=q_add){
                int query_rank = (int) (q * queryN);

                double merge_v = longToResult(merge_sst_worker.findMinValueWithRank(query_rank));
                int merge_delta_rank = getDeltaRank(query_a, queryN, merge_v, query_rank);
                double merge_relative_err = 1.0*merge_delta_rank/(queryN);
                err_merge_sst+=Math.abs(merge_relative_err)/(q_count*TEST_CASE);

//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("sig="+sig+"\t\t" + queryN + "\t" + "\t\t\t\t\t\t" + "\t" + err_merge_sst+"\t"+1.0*mergedSketchN/TEST_CASE);
        //        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_mergeBuf-err_full)+"\t"+queryN*(err_merge-err_full));
        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t\t\t\t\t" + "\t" + err_merge_sst+"\t"+1.0*mergedSketchN/TEST_CASE));
//        time_result.set(RESULT_LINE, time_result.get(RESULT_LINE).concat("\t\t\t\t\t\t" + merge_time+(TEST_FULL?("\t"+full_time):"")));

    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START=new Date().getTime();
        IntervalEvaluatingLatencySSTKLL main;

        System.out.println("UnSeqData. LSM-KLL. interval query" + "\n");
        for (int startType=1,endType=1,dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            main = new IntervalEvaluatingLatencySSTKLL();
            main.prepareA(dataType);
            for (double mu : muList)
                for (double sig : sigList) {
                    main.prepareDisorder(mu, sig);
                    for (int queryN : new int[]{40000000})
                        for (int query_mem : new int[]{1024 * 256})
                            for (int chunk_seri : new int[]{/*128,256,512,*/4096/*,128,256,512,1024/*,2048,4096,8192*/})
                                for (int SketchSizeRatio : new int[]{/*1,2,4,8,16,32*/4}) {
                                    if (dataType == startType) {
                                        err_result.add("sig:\t"+sig+"\tN:" + queryN + ", " + "T_s:" + SketchSizeRatio + ", " + "M:\t" + query_mem / 1024 + "\t");
                                        time_result.add("N:" + queryN + ", " + "T_s:" + SketchSizeRatio + ", " + "M:\t" + query_mem / 1024 + "\t");
                                    }
                                    main.testMergeError(queryN, query_mem, chunk_seri, sig, SketchSizeRatio);
                                    main.RESULT_LINE++;
                                }
                }
        }
        System.out.println("LSM KLL\nTEST_CASE=" + TEST_CASE);
        System.out.println("\nError rate:");
        for (String s : err_result)
            System.out.println(s);
        System.out.println("\t\tALL_TIME:"+(new Date().getTime()-START));
    }
}
