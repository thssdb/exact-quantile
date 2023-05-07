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
public class MainForKLLDivideBound {
    int dataType=-233;
    static double Pr = 0.5;
    static int pagePerChunkL=16,pagePerChunkR=16;
    static int startType=2,endType=2;
    static int chunkN/* = 4096*16 */, chunkNum/* = N/chunkN*/, pageN = 4096;
    static int N = pageN * (4096), pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 1; // CHECK IT
    static double[] a;
    static KLLSketchForQuantile[] KLLWorkerByChunkDivide;
    static long[] workerMinT, workerMaxT,chunkMinV,chunkMaxV;
    static int compaction_level;
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
//                    random.nextDouble();
//                    i-0.5*Math.floor(i/pageN)*pageN;
                    longToResult(i/32768);
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

    public void prepareWorker(int maxSeriesByte) {
        compaction_level = 0;
        for(int i=0;i<=10;i++)if((maxSeriesByte<<i)>=pageN*8){
            compaction_level=i;
            break;
        }
        KLLWorkerByChunkDivide= new KLLSketchForQuantile[pageNum];
        workerMinT = new long[pageNum];
        workerMaxT = new long[pageNum];
        chunkMinV = new long[chunkNum];
        chunkMaxV = new long[chunkNum];
        int enoughMemByte = pageN * 9;
        for (int i = 0; i < pageNum; i++) {
            workerMinT[i] = Long.MAX_VALUE;
            workerMaxT[i] = Long.MIN_VALUE;
            for (int j = 0; j < pageN; j++) {
                int index = (i * pageN + j);
                workerMinT[i] = Math.min(workerMinT[i], index);
                workerMaxT[i] = Math.max(workerMaxT[i], index);
            }
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
//            System.out.print("\t\tPageMinMaxValue:");
//            for(int j=0;j<pageNumInChunk;j++)System.out.print("\t"+longToResult(pageMinV.getLong(j))+"..."+longToResult(pageMaxV.getLong(j))+"\t");
//            System.out.println();
//            System.out.println("\t\tSortedPageValue:");
//            for(int j=0;j<pageNumInChunk;j++){
//                System.out.print("\t\t\t\t");
//                double[] tmp = new double[pageN];
//                for(int k=j*pageN;k<j*pageN+pageN;k++){
//                    tmp[k-j*pageN]=longToResult(v.getLong(k));
//                }
//                Arrays.sort(tmp);
//                for(int k=0;k<pageN;k++)System.out.print(tmp[k]+" ");
//                System.out.println();
//            }
            KLLSketchDividedForBound chunkSketch = new KLLSketchDividedForBound(compaction_level,v);
//            for(int j=i*chunkN/pageN;j<i*chunkN/pageN+chunkN/pageN;j++)
//                KLLWorkerByChunkDivide[j] = new LongKLLSketchForUnit(chunkSketch,workerMaxT[j],pageN);
            List<LongKLLSketchForUnit> pageSketch = chunkSketch.divideMemSketchByItemValue(pageStartIndex,pageMinV,pageMaxV);
            chunkMinV[i] = chunkSketch.minV;
            chunkMaxV[i] = chunkSketch.maxV;
//            System.out.println("\t\t\tchunkExtremeValue:\t"+longToResult(chunkMinV[i])+"\t\t"+longToResult(chunkMaxV[i]));
            for(int j=0;j<pageNumInChunk;j++) {
                KLLWorkerByChunkDivide[i * pageNumInChunk + j] = pageSketch.get(j);
//                if(i<=5) {
//                    pageSketch.get(j).show();
//                    pageSketch.get(j).showNum();
//                }
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
        double err_full = 0, err_merge_page = 0, err_merge_divide=0;
        double prOutOfBound = 0;
        double[] query_a = new double[queryN];
        assert queryN%chunkN==0;

        for (int T = 0; T < TEST_CASE; T++) {
            System.out.println("\tTEST_CASE:"+T);
            prepareWorker(maxSeriesByte);
//            System.out.println("????");
//            KLLWorkerByPage[0].show();
//            KLLWorkerByChunkDivide[0].show();
//            int rkA=KLLWorkerByPage[0].getApproxRank(dataToLong(-0.0029907192953519));
//            int rkB=KLLWorkerByChunkDivide[0].getApproxRank(dataToLong(-0.0029907192953519));
//
//            System.out.println("!!!!!!!!!!!!!!!!!!!!!!rkAB:"+rkA+"  "+rkB);
            for (int L = 0; L + queryN <= N; L += queryN) {
                int R = L+queryN;
//            System.out.println("\t\t\t"+posL+"\t"+posR);

                merge_page_time -= new Date().getTime();
//                HeapLongStrictKLLSketch merge_page_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
//                for (int i = L / pageN; i < R / pageN; i++)
//                    merge_page_worker.mergeWithTempSpace(Collections.singletonList(KLLWorkerByPage[i]));

                FastKLLSketchLazyForBound merge_divide_worker = new FastKLLSketchLazyForBound(maxMemoryByte);
                for (int i = L / pageN; i < R / pageN; i++)
                    merge_divide_worker.mergeWithTempSpace(Collections.singletonList(KLLWorkerByChunkDivide[i]));
                for(int chunkID=L/chunkN;chunkID<R/chunkN;chunkID++){
                    merge_divide_worker.addRecord(true,chunkMinV[chunkID],chunkMaxV[chunkID],compaction_level);
                }

                if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
                Arrays.sort(query_a);

                double q_start = 0.01, q_end = 0.99, q_add = 0.01, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
                double ratioPerQuery = 1.0/(q_count * TEST_CASE*N/queryN);
                for (double q = q_start; q < q_end + 1e-10; q += q_add) {

                    int query_rank = (int) (q * queryN);

                    double merge_divide_v = longToResult(merge_divide_worker.findMinValueWithRank(query_rank));
                    int merge_divide_delta_rank = getValueActualRank(query_a, queryN, merge_divide_v) - query_rank;
                    double merge_divide_relative_err = 1.0 * merge_divide_delta_rank / (queryN);
                    err_merge_divide += Math.abs(merge_divide_relative_err) *ratioPerQuery;

                    double bound = merge_divide_worker.queryBound(dataToLong(merge_divide_v),Pr)/queryN;

                    if(merge_divide_relative_err>bound)prOutOfBound+=ratioPerQuery;

                    System.out.println("Q="+q+"\tapprox=\t"+merge_divide_v+"\t\tbound\t"+bound+"\t\tactErr\t"+merge_divide_relative_err);

//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
//                    System.out.println("\t\t\trk:"+query_rank+"\t\t"+Math.abs(merge_divide_relative_err));
                }
//                System.out.println("\n");
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("TEST_CASE="+TEST_CASE+"\t\t" + queryN +"\tPagesInChunk:"+chunkN/pageN+ "\t" + err_merge_page + "\t"+err_merge_divide+"\t\tprOutOfBound\t"+prOutOfBound+"\t\tPr\t"+Pr);
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        MainForKLLDivideBound main;

        System.out.println("TEST for Correct Sketch Divide" + "\n");
        main = new MainForKLLDivideBound();
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            main.RESULT_LINE=0;
            main.prepareA(dataType);
            for (int PagePerChunk=pagePerChunkL;PagePerChunk<=pagePerChunkR;PagePerChunk*=2) {
                chunkN=pageN*PagePerChunk;
                chunkNum = N/chunkN;
                if(dataType==startType)
                    main.err_result.add("Pages in Chunk:\t\t" + PagePerChunk + "\t");
                for (int i : new int[]{/*pageN,chunkN,*/N})
                    for (int query_mem : new int[]{/*1024*16,1024*32,1024*64,*/1024 * 128/*,1024*256*//*,1024*512,1024*1024*/})
                        for (int page_seri : new int[]{/*128,256,512,*/pageN*8/32/*,128,256,512,1024/*,2048,4096,8192*/})
                            main.testError(i, query_mem, page_seri);
                main.RESULT_LINE++;
            }
            System.out.println("byPage & byChunkDivide\nTEST_CASE=" + TEST_CASE);
            System.out.println("\nError rate:");
        }
        for (String s : main.err_result)
            System.out.println(s);
    }
}
