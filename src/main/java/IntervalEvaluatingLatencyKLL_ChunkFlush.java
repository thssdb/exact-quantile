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

public class IntervalEvaluatingLatencyKLL_ChunkFlush {
    int dataType = 1;
    static int N = 4096 * 5000; // CHECK IT
    static int chunkN = 4096*16 , chunkNum = N/chunkN, pageN = 4096, pageNum = N / pageN;
    public static int TEST_CASE = 1024; // CHECK IT
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
    static double mu = 2.0, sig = 1.5;

    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] = Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
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
    }

    public void prepareWorker(int maxSeriesByte) {
        int compactionLevel = 0;
        for(int i=0;i<=10;i++)if((maxSeriesByte<<i)>=pageN*8){
            compactionLevel=i;
            break;
        }
        KLLWorkerByPage = new KLLSketchForQuantile[pageNum];
        KLLWorkerByChunkDivide= new KLLSketchForQuantile[pageNum];
        workerMinT = new long[pageNum];
        workerMaxT = new long[pageNum];
        int enoughMemByte = pageN * 9;
        for (int i = 0; i < pageNum; i++) {
            LongKLLSketch worker = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);

            workerMinT[i] = Long.MAX_VALUE;
            workerMaxT[i] = Long.MIN_VALUE;
            for (int j = 0; j < pageN; j++) {
                int index = cc.getInt(i * pageN + j);
                worker.update(dataToLong(a[index]));
                workerMinT[i] = Math.min(workerMinT[i], index);
                workerMaxT[i] = Math.max(workerMaxT[i], index);
            }
            worker.compactBeforeSerialization();
            KLLWorkerByPage[i] = worker;
//            if(i==0)worker.show();
        }
        for(int i=0;i<chunkNum;i++){
            ObjectArrayList<IntLongPair> TVList = new ObjectArrayList<>();
            for(int j=i*chunkN;j<i*chunkN+chunkN;j++){
                int index = cc.getInt(j);
                TVList.add(PrimitiveTuples.pair(index,dataToLong(a[index])));
            }
            TVList.sort(Comparator.comparingInt(IntLongPair::getOne));
            LongArrayList t=new LongArrayList(chunkN),v=new LongArrayList((chunkN));
            for(int j=0;j<chunkN;j++){
                t.add(TVList.get(j).getOne());
                v.add(TVList.get(j).getTwo());
            }

            int pageNumInChunk =chunkN/pageN;
            IntArrayList pageStartIndex = new IntArrayList(pageNumInChunk);
            LongArrayList pageMinV = new LongArrayList(pageNumInChunk),pageMaxV = new LongArrayList(pageNumInChunk);
            for(int j=0;j<pageNumInChunk;j++) {
                pageStartIndex.add(j * pageN);
                long minV=Long.MAX_VALUE,maxV=Long.MIN_VALUE;
                for(int k=i*chunkN+j*pageN;k<i*chunkN+j*pageN+pageN;k++){
                    int index = k;
                    minV=Math.min(minV,dataToLong(a[index]));
                    maxV=Math.max(maxV,dataToLong(a[index]));
                }
                pageMinV.add(minV);
                pageMaxV.add(maxV);
            }
            KLLSketchForDivide chunkSketch = new KLLSketchForDivide(compactionLevel,t,v);
//            for(int j=i*chunkN/pageN;j<i*chunkN/pageN+chunkN/pageN;j++)
//                KLLWorkerByChunkDivide[j] = new LongKLLSketchForUnit(chunkSketch,workerMaxT[j],pageN);
            List<LongKLLSketchForUnit> pageSketch = chunkSketch.divideMemSketchByItemValue(pageStartIndex,pageMinV,pageMaxV);
            for(int j=0;j<pageNumInChunk;j++) {
                KLLWorkerByChunkDivide[i * pageNumInChunk + j] = pageSketch.get(j);
//                pageSketch.get(j).showNum();
            }
        }
//        for(int i=0;i<100;i++)System.out.println(workerMinT[i]+"---"+workerMaxT[i]);
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
        return !inInterval(x, y, L, R) && !(y < L || x > R);
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
        IntervalEvaluatingLatencyKLL_ChunkFlush main;

        System.out.println("interval query" + "\n");
        for (int dataType = 1; dataType <= 1; dataType++) { // CHECK IT
            main = new IntervalEvaluatingLatencyKLL_ChunkFlush();
            main.prepareA(dataType);
            for (double mu = 2.0, sig = 0.0; sig < 0.1; sig += 0.5) {
                IntervalEvaluatingLatencyKLL_ChunkFlush.mu = mu;
                IntervalEvaluatingLatencyKLL_ChunkFlush.sig = sig;
                main.prepareDisorder();
                main.err_result.add("mu:" + mu + ", " + "sig:" + sig + "\t");
                main.time_result.add("mu:" + mu + ", " + "sig:" + sig + "\t");
                for (int i : new int[]{100000,1000000,10000000})
                    for (int query_mem : new int[]{/*1024*16,1024*32,1024*64,*/1024 * 128/*,1024*256*//*,1024*512,1024*1024*/})
                        for (int page_seri : new int[]{/*128,256,512,*/1024/*,128,256,512,1024/*,2048,4096,8192*/})
                            main.testMergeError(i, query_mem, page_seri);
                main.RESULT_LINE++;
            }
            System.out.println("LSM-KLL & KLL\nTEST_CASE=" + TEST_CASE);
            System.out.println("\nError rate:");
            for (String s : main.err_result)
                System.out.println(s);
            System.out.println("\nQuery Time:");
            for (String s : main.time_result)
                System.out.println(s);
        }
    }
}
