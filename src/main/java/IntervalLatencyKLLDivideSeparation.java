import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.eclipse.collections.api.tuple.primitive.IntLongPair;
import org.eclipse.collections.impl.tuple.primitive.PrimitiveTuples;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class IntervalLatencyKLLDivideSeparation {
    int dataType = 0;
    static int N = /*50000000/4096*4096*/16*12,sketchCompression=4; // CHECK IT
    static int chunkN = /*4096*32*/16*2, chunkNum = N/chunkN, pageN = 16, pageNum = N / pageN;
    public static int TEST_CASE = 1; // CHECK IT
//    static int timedFlushUnSeq = ??;
    boolean TEST_FULL = false; // CHECK IT
    static double[] a;
    static KLLSketchForQuantile[] KLLWorkerByPage,KLLWorkerByChunkDivide;
    static long[] workerMinT, workerMaxT;
     ArrayList<String> err_result = new ArrayList<>();
     ArrayList<String> time_result = new ArrayList<>();
    boolean show_time = false, show_err = true;
     int RESULT_LINE = 0;
    Random random = new Random(233);
    static long[] aa;
    static long[] bb;
    static IntArrayList cc;
    static double mu = 8.0, sig = 4.0;


    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] =
                    //Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
//                    i;
            longToResult(i);
            return;
        }
        BufferedReader reader = null;
        if (dataType == 1)
            reader = new BufferedReader(new FileReader(new File("1_bitcoin.csv")));
        if (dataType == 2)
            reader = new BufferedReader(new FileReader(new File("2_physiological_stress.txt")));
        if (dataType == 3)
            reader = new BufferedReader(new FileReader(new File("4_taxipredition8M.txt")));
        if (dataType == 4)
            reader = new BufferedReader(new FileReader(new File("5_wh.csv")));
        assert reader != null;
        reader.readLine(); // ignore first line.
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if (cntN == N) break;
        }
    }

    public void prepareDisorder() {
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
        for(int i=0;i<N;i++){
            System.out.print("\t\t"+cc.getInt(i)+" "+bb[cc.getInt(i)]+" "+a[cc.getInt(i)]);
            if((i+1)%pageN==0)System.out.println();
        }

    }


    public void crossCompaction(int compactionLevel, ObjectArrayList<IntLongPair> unSeqTVList,ObjectArrayList<LongKLLSketchForUnit> pageSketch, ObjectArrayList<ObjectArrayList<IntLongPair>> pageData, LongArrayList pageMinT,LongArrayList pageMaxT){
        unSeqTVList.sort(Comparator.comparingInt(IntLongPair::getOne));
        long minT=unSeqTVList.get(0).getOne(),maxT=unSeqTVList.get(unSeqTVList.size()-1).getOne();
        int pageNum = pageSketch.size();
        ObjectArrayList<IntLongPair> allTVList = new ObjectArrayList<>(unSeqTVList);
        boolean[] overlapped=new boolean[pageNum];
        System.out.println("\tcrossCompaction:\tunSeqTVSize:"+unSeqTVList.size()+"\tunSeqTime:"+minT+"..."+maxT+"\t\tflushedPageNum:"+pageSketch.size());
        int overlappedPageNum=0;
        for(int i=0;i<pageNum;i++)if(overlapInterval(pageMinT.getLong(i),pageMaxT.getLong(i),minT,maxT)){
            allTVList.addAll(pageData.get(i));
            overlapped[i]=true;
            overlappedPageNum++;
//            System.out.println("\t\t\tOverlapped pageT"+pageMinT.getLong(i)+"..."+pageMaxT.getLong(i));
        }
        System.out.println("\t\t\tOverlappedPageNum:"+overlappedPageNum);
        for(int i=pageNum-1;i>=0;i--)if(overlapped[i]) {
            pageSketch.remove(i);
            pageData.remove(i);
            pageMinT.removeLong(i);
            pageMaxT.removeLong(i);
        }
//        System.out.println("|DEBUG\t\t\t allTVSize"+allTVList.size());
        getMemSketchAndDivide(compactionLevel,allTVList,pageSketch,pageData,pageMinT,pageMaxT);
    }

    public void finalCompaction(int compactionLevel, ObjectArrayList<IntLongPair> unSeqTVList,ObjectArrayList<IntLongPair> seqTVList,ObjectArrayList<LongKLLSketchForUnit> pageSketch, ObjectArrayList<ObjectArrayList<IntLongPair>> pageData, LongArrayList pageMinT,LongArrayList pageMaxT){
        unSeqTVList.sort(Comparator.comparingInt(IntLongPair::getOne));
        long minT=unSeqTVList.get(0).getOne(),maxT=unSeqTVList.get(unSeqTVList.size()-1).getOne();
        int pageNum = pageSketch.size();
        ObjectArrayList<IntLongPair> allTVList = new ObjectArrayList<>(unSeqTVList);
        boolean[] overlapped=new boolean[pageNum];
        System.out.println("\tfinalCrossCompaction:\tunSeqTVSize:"+unSeqTVList.size()+"\tunSeqTime:"+minT+"..."+maxT+"\t\tflushedPageNum:"+pageSketch.size()+"\t\tseqTVSize:"+seqTVList.size());
        for(int i=0;i<pageNum;i++)if(overlapInterval(pageMinT.getLong(i),pageMaxT.getLong(i),minT,maxT)){
            allTVList.addAll(pageData.get(i));
            overlapped[i]=true;
        }
        for(int i=pageNum-1;i>=0;i--)if(overlapped[i]) {
            pageSketch.remove(i);
            pageData.remove(i);
            pageMinT.removeLong(i);
            pageMaxT.removeLong(i);
        }
        allTVList.addAll(seqTVList);
        getMemSketchAndDivide(compactionLevel,allTVList,pageSketch,pageData,pageMinT,pageMaxT);
    }

    private void getMemSketchAndDivide(int compactionLevel, ObjectArrayList<IntLongPair> TVList,ObjectArrayList<LongKLLSketchForUnit> pageSketch, ObjectArrayList<ObjectArrayList<IntLongPair>> pageData, LongArrayList pageMinT,LongArrayList pageMaxT){
        int memN = TVList.size();
        TVList.sort(Comparator.comparingInt(IntLongPair::getOne));
        IntArrayList t=new IntArrayList(memN);
        LongArrayList v=new LongArrayList(memN);
        for(int j=0;j<memN;j++){
            t.add(TVList.get(j).getOne());
            v.add(TVList.get(j).getTwo());
        }
        TVList.clear();
        TVList.trim(10);
        int pageNumInChunk =memN/pageN;
        IntArrayList pageStartIndex = new IntArrayList(pageNumInChunk);
        LongArrayList pageMinV = new LongArrayList(pageNumInChunk),pageMaxV = new LongArrayList(pageNumInChunk);
        for(int j=0;j<pageNumInChunk;j++) {
            pageStartIndex.add(j * pageN);
            long minV=Long.MAX_VALUE,maxV=Long.MIN_VALUE;
            for(int k=j*pageN;k<j*pageN+pageN;k++){
                minV=Math.min(minV,v.getLong(k));
                maxV=Math.max(maxV,v.getLong(k));
            }
            pageMinV.add(minV);
            pageMaxV.add(maxV);
        }
        KLLSketchForDivide memSketch = new KLLSketchForDivide(compactionLevel,v);
        List<LongKLLSketchForUnit> pageSketchInChunk = memSketch.divideMemSketchByItemValue(pageStartIndex,pageMinV,pageMaxV);
//        ObjectArrayList<Pair<LongKLLSketchForUnit,ObjectArrayList<IntLongPair>>> result= new ObjectArrayList<>();
        for(int j=0;j<pageNumInChunk;j++){
            ObjectArrayList<IntLongPair> pageTV= new ObjectArrayList<>();
            for(int i=j*pageN;i<j*pageN+pageN;i++)
                pageTV.add(PrimitiveTuples.pair((int)t.getInt(i),v.getLong(i)));
            pageData.add(pageTV);
            pageSketch.add(pageSketchInChunk.get(j));
            pageMinT.add(t.getInt(j*pageN));
            pageMaxT.add(t.getInt(j*pageN+pageN-1));
        }
    }

    public void prepareWorker(int maxSeriesByte) {
        int compactionLevel = 0;
        for(int i=0;i<=10;i++)if((maxSeriesByte<<i)>=pageN*8){
            compactionLevel=i;
            break;
        }

        ObjectArrayList<LongKLLSketchForUnit> pageSketch=new ObjectArrayList<>();
        ObjectArrayList<ObjectArrayList<IntLongPair>> pageData=new ObjectArrayList<>();
        LongArrayList pageMinT=new LongArrayList(),pageMaxT=new LongArrayList();
        ObjectArrayList<IntLongPair> seqTVList = new ObjectArrayList<>();
        ObjectArrayList<IntLongPair> unSeqTVList = new ObjectArrayList<>();
        long lastFlushSeqEndTime = Long.MIN_VALUE;
        for(int i=0;i<N;i++){
            int index = cc.getInt(i);
            int cntT = index;
            long cntV = dataToLong(a[index]);
            if(cntT<=lastFlushSeqEndTime){
                unSeqTVList.add(PrimitiveTuples.pair(cntT,cntV));
                if(unSeqTVList.size()==chunkN){//不会影响最新的seq
//                    for(int k=0;k<pageSketch.size();k++)System.out.println("\t\t\t\tPageT:"+pageMinT.getLong(k)+"..."+pageMaxT.getLong(k));
                    crossCompaction(compactionLevel,unSeqTVList,pageSketch,pageData,pageMinT,pageMaxT);
                    unSeqTVList.clear();
//                    break;
                }
            }else {
                seqTVList.add(PrimitiveTuples.pair(cntT,cntV));
                if(seqTVList.size()==chunkN){
                    getMemSketchAndDivide(compactionLevel,seqTVList,pageSketch,pageData,pageMinT,pageMaxT);
                    lastFlushSeqEndTime = pageMaxT.getLong(pageMaxT.size()-1);
                    seqTVList.clear();
                }
            }
        }
        if(unSeqTVList.size()>0){ // then seqTVList.size()>0, since N=pageN*?
            finalCompaction(compactionLevel,unSeqTVList,seqTVList,pageSketch,pageData,pageMinT,pageMaxT);
        }

        KLLWorkerByChunkDivide= new KLLSketchForQuantile[pageNum];
        IntArrayList pageIndex = new IntArrayList();
        for(int i=0;i<pageNum;i++)pageIndex.add(i);
        pageIndex.sort((x,y)->Long.compare(pageMinT.getLong(x),pageMinT.getLong(y)));
        for(int i=0;i<pageNum;i++){
            KLLWorkerByChunkDivide[i]=pageSketch.get(pageIndex.getInt(i));
            KLLWorkerByChunkDivide[i].show();
            KLLWorkerByChunkDivide[i].showNum();
        }
//        for(int i=0;i<pageSketch.size();i++)System.out.println("\t\t\t\tPageT:"+pageMinT.getLong(i)+"..."+pageMaxT.getLong(i)+"\tpageSize:"+pageData.get(i).size());
    }

    public int getValueActualRank(double[] sortedA, int queryN, double v) { // number of elements <= v
        int L = 0, R = queryN - 1;
        while (L < R) {
            int mid = (L + R + 1) >>> 1;
            if (v < sortedA[mid]) R = mid - 1;
            else L = mid;
        }
        return L;
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
        return /*!inInterval(x, y, L, R) && */!(y < L || x > R);
    }


    public void testMergeError(int queryN, int maxMemoryByte, int maxSeriesByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time = 0, merge_time = 0;
        double err_full = 0, err_merge = 0;
        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        double  avgAvailable=0;
        for (int T = 0; T < TEST_CASE; T++) {
            prepareWorker(maxSeriesByte);
            int L = LL[T], R = RR[T];
//            System.out.println("\t\t\t"+posL+"\t"+posR);

            merge_time -= new Date().getTime();
            int buf_kll_num = 1;
            List<KLLSketchForQuantile> buf_kll_list = new ArrayList<>(buf_kll_num);
            HeapLongStrictKLLSketch merge_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            int available = 0;
            for (int i = 0; i < pageNum; i++)
                if (inInterval(workerMinT[i], workerMaxT[i], L, R - 1)) { // [L,R)
                    available++;
                    if (buf_kll_list.size() < buf_kll_num)
                        buf_kll_list.add(KLLWorkerByPage[i]);
                    else {
                        merge_worker.mergeWithTempSpace(buf_kll_list);
                        buf_kll_list.clear();
                        buf_kll_list.add(KLLWorkerByPage[i]);
                    }
                } else if (overlapInterval(workerMinT[i], workerMaxT[i], L, R - 1)) {
                    for (int j = i * pageN; j < (i + 1) * pageN; j++) {
                        int index = cc.getInt(j);
                        if (inInterval(index, L, R - 1))
                            merge_worker.update(dataToLong(a[index]));
                    }
                }
            avgAvailable+=1.0*available*pageN/queryN/TEST_CASE;
//            System.out.println("\t\tavailable:"+1.0*available*pageN/queryN);

            full_time -= new Date().getTime();
            HeapLongStrictKLLSketch full_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            if (TEST_FULL) {
                for (int i = L; i < R; i++)
                    full_worker.update(dataToLong(a[i]));
            }
            full_time += new Date().getTime();

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);

            double q_start = 0.01, q_end = 0.99, q_add = 0.0025, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
                int query_rank = (int) (q * queryN);

                double merge_v = longToResult(merge_worker.findMinValueWithRank(query_rank));
                int merge_delta_rank = getValueActualRank(query_a, queryN, merge_v) - query_rank;
                double merge_relative_err = 1.0 * merge_delta_rank / (queryN);
                err_merge += Math.abs(merge_relative_err) / (q_count * TEST_CASE);

                if (TEST_FULL) {
                    double full_v = longToResult(full_worker.findMinValueWithRank(query_rank));
                    int full_delta_rank = getValueActualRank(query_a, queryN, full_v) - query_rank;
                    double full_relative_err = 1.0 * full_delta_rank / (queryN);
                    err_full += Math.abs(full_relative_err) / (q_count * TEST_CASE);
                }

//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("\t\t" + queryN + "\t" + err_merge + (TEST_FULL ? ("\t" + err_full) : "")+"\t\tavgAvailable:"+avgAvailable);
        //        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_mergeBuf-err_full)+"\t"+queryN*(err_merge-err_full));
        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t\t\t\t\t" + err_merge + (TEST_FULL ? ("\t" + err_full) : "")));
//        time_result.set(RESULT_LINE, time_result.get(RESULT_LINE).concat("\t\t\t\t\t\t" + merge_time+(TEST_FULL?("\t"+full_time):"")));
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        IntervalLatencyKLLDivideSeparation main;

        System.out.println("interval query" + "\n");
        for (int dataType = 0; dataType <= 0; dataType++) { // CHECK IT
            main = new IntervalLatencyKLLDivideSeparation();
            main.prepareA(dataType);
            main.prepareDisorder();
            main.prepareWorker(pageN*8/sketchCompression);
//            for (double mu = 2.0, sig = 0.0; sig < 0.1; sig += 0.5) {
//                IntervalLatencyKLLDivideSeparation.mu = mu;
//                IntervalLatencyKLLDivideSeparation.sig = sig;
//                main.prepareDisorder();
//                main.err_result.add("mu:" + mu + ", " + "sig:" + sig + "\t");
//                main.time_result.add("mu:" + mu + ", " + "sig:" + sig + "\t");
//                for (int i : new int[]{100000,1000000,10000000})
//                    for (int query_mem : new int[]{/*1024*16,1024*32,1024*64,*/1024 * 128/*,1024*256*//*,1024*512,1024*1024*/})
//                        for (int page_seri : new int[]{/*128,256,512,*/1024/*,128,256,512,1024/*,2048,4096,8192*/})
//                            main.testMergeError(i, query_mem, page_seri);
//                main.RESULT_LINE++;
//            }
//            System.out.println("LSM-KLL & KLL\nTEST_CASE=" + TEST_CASE);
//            System.out.println("\nError rate:");
//            for (String s : main.err_result)
//                System.out.println(s);
//            System.out.println("\nQuery Time:");
//            for (String s : main.time_result)
//                System.out.println(s);
        }
    }
}
