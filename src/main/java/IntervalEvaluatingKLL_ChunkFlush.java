import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
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

/** Compare sketchByPage & sketchByChunkDivide.
 * Different ChunkSize.
*/
public class IntervalEvaluatingKLL_ChunkFlush {
    int dataType=-233;
    static int pagePerChunkL=32,pagePerChunkR=32;
    static int chunkN/* = 4096*16 */, chunkNum/* = N/chunkN*/, pageN = 4096;
    static int N = 50000000/4096*4096, pageNum = N / pageN; // CHECK IT  Can't query 5E7.
    public static int TEST_CASE = 512; // CHECK IT
    static double[] a;
    static KLLSketchForQuantile[] KLLWorkerByPage,KLLWorkerByChunkDivide;
    static long[] workerMinT, workerMaxT;
    ArrayList<String> err_result = new ArrayList<>();
    ArrayList<String> time_result = new ArrayList<>();
    boolean show_time = false, show_err = true;
    int RESULT_LINE = 0;
    Random random = new Random(233);

    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] = //longToResult(i&4095);
                        //Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
                    random.nextDouble();
            return;
        }
        BufferedReader reader = null;
        if (dataType == 1)
            reader = new BufferedReader(new FileReader(new File("1_bitcoin.csv")));
        if (dataType == 2)
            reader = new BufferedReader(new FileReader(new File("SpacecraftThruster.txt")));
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
                int index = (i * pageN + j);
                worker.update(dataToLong(a[index]));
                workerMinT[i] = Math.min(workerMinT[i], index);
                workerMaxT[i] = Math.max(workerMaxT[i], index);
            }
            worker.compactBeforeSerialization();
            KLLWorkerByPage[i] = worker;
//            worker.showNum();
        }
        for(int i=0;i<chunkNum;i++){
            ObjectArrayList<IntLongPair> TVList = new ObjectArrayList<>();
            for(int j=i*chunkN;j<i*chunkN+chunkN;j++){
                int index = j;
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
                for(int k=j*pageN;k<j*pageN+pageN;k++){
                    minV=Math.min(minV,v.getLong(k));
                    maxV=Math.max(maxV,v.getLong(k));
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
//            if(i==0)KLLWorkerByChunkDivide[0].showNum();
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


    public void testError(int queryN, int maxMemoryByte, int maxSeriesByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time = 0, merge_page_time = 0;
        double err_full = 0, err_merge_page = 0, err_merge_divide = 0;
        double[] query_a = new double[queryN];
        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        DoubleArrayList errPage=new DoubleArrayList(),errDivide=new DoubleArrayList(),errFull=new DoubleArrayList();

        for (int T = 0; T < TEST_CASE; T++) {
            if((T&1)==0)//CHECK IT!
                prepareWorker(maxSeriesByte);
            int L = LL[T], R = RR[T];
            int pageL = (L+pageN-1)/pageN, pageR = R/pageN;
            int posL = pageL*pageN, posR = pageR*pageN;
//            System.out.println("\t\t\t"+posL+"\t"+posR);

            merge_page_time -= new Date().getTime();
            HeapLongStrictKLLSketch merge_page_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            for(int i=L;i<Math.min(R,posL);i++)
                merge_page_worker.update(dataToLong(a[i]));
            for (int i=pageL;i<pageR; i++)
                merge_page_worker.mergeWithTempSpace(Collections.singletonList(KLLWorkerByPage[i]));
            for(int i=Math.max(L,posR);i<R;i++)
                merge_page_worker.update(dataToLong(a[i]));

            HeapLongStrictKLLSketch merge_divide_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            for(int i=L;i<Math.min(R,posL);i++)
                merge_divide_worker.update(dataToLong(a[i]));
            for (int i=pageL;i<pageR; i++)
                merge_divide_worker.mergeWithTempSpace(Collections.singletonList(KLLWorkerByChunkDivide[i]));
            for(int i=Math.max(L,posR);i<R;i++)
                merge_divide_worker.update(dataToLong(a[i]));

            HeapLongStrictKLLSketch full_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            for(int i=L;i<R;i++)
                full_worker.update(dataToLong(a[i]));

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);

            double q_start = 0.01, q_end = 0.99, q_add = 0.001, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
            double ratioPerQuery = 1.0 / (q_count * TEST_CASE);
            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
                int query_rank = (int) (q * queryN);

                double merge_page_v = longToResult(merge_page_worker.findMinValueWithRank(query_rank));
                int merge_page_delta_rank = getValueActualRank(query_a, queryN, merge_page_v) - query_rank;
                double merge_page_relative_err = 1.0 * merge_page_delta_rank / (queryN);
                err_merge_page += Math.abs(merge_page_relative_err) * ratioPerQuery;

                double merge_divide_v = longToResult(merge_divide_worker.findMinValueWithRank(query_rank));
                int merge_divide_delta_rank = getValueActualRank(query_a, queryN, merge_divide_v) - query_rank;
                double merge_divide_relative_err = 1.0 * merge_divide_delta_rank / (queryN);
                err_merge_divide += Math.abs(merge_divide_relative_err) *ratioPerQuery;

                double full_v = longToResult(full_worker.findMinValueWithRank(query_rank));
                int full_delta_rank = getValueActualRank(query_a, queryN, full_v) - query_rank;
                double full_relative_err = 1.0 * full_delta_rank / (queryN);
                err_full += Math.abs(full_relative_err) *ratioPerQuery;

                errPage.add(Math.abs(merge_page_relative_err));
                errDivide.add(Math.abs(merge_divide_relative_err));
                errFull.add(Math.abs(full_relative_err));
//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
//                    System.out.println("\t\t\trk:"+query_rank+"\t\t"+Math.abs(merge_divide_relative_err));
            }

        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("\t\t" + queryN + "\t" + err_merge_page + "\t" + err_merge_divide + "\t" + err_full);
        //        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_mergeBuf-err_full)+"\t"+queryN*(err_merge-err_full));
        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t\t\t\t\t" + err_merge_page + "\t" + err_merge_divide + "\t" + err_full));

        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t|\t" + calcSTDDEV(errPage) + "\t" + calcSTDDEV(errDivide) + "\t" + calcSTDDEV(errFull)+"\t|\t"));
    }

    private double calcSTDDEV(DoubleArrayList err){
        double avg = 0;
        for(double x:err)avg+=x;
        avg/=err.size();
        double STDDEV=0;

        for(double x:err)STDDEV+=(x-avg)*(x-avg);
        STDDEV/=err.size();
        STDDEV = Math.sqrt(STDDEV);
        return STDDEV;
    }


    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START=new Date().getTime();
        IntervalEvaluatingKLL_ChunkFlush main;

        System.out.println("range query: byPage & byChunkDivide" + "\n");
        for (int dataType = 1; dataType <= 2; dataType++) { // CHECK IT
            main = new IntervalEvaluatingKLL_ChunkFlush();
            main.prepareA(dataType);
            for (int PagePerChunk=pagePerChunkL;PagePerChunk<=pagePerChunkR;PagePerChunk*=2)
                for (int queryN : new int[]{100000,400000,1000000,4000000,10000000,20000000,40000000}){
                chunkN=pageN*PagePerChunk;
                chunkNum = N/chunkN;
                main.err_result.add("Pages_in_Chunk:\t" + PagePerChunk + "\tqueryN:\t"+queryN);
                    for (int query_mem : new int[]{1024 * 128})
                        for (int page_seri : new int[]{1024})
                            main.testError(queryN, query_mem, page_seri);
                main.RESULT_LINE++;
            }
            System.out.println("byPage & byChunkDivide\nTEST_CASE=" + TEST_CASE);
            System.out.println("\nError rate:           err_mergePage    err_divideChunk   err_NoPreComputation |   STD DEV * 3");
            for (String s : main.err_result)
                System.out.println(s);
            System.out.println("\n------------------------------------------------------------------\n\n");
        }
        System.out.println("\t\tALL_TIME:"+(new Date().getTime()-START));
    }
}
