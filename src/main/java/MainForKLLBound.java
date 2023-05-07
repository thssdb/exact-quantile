import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

/** Compare sketchByPage & sketchByChunkDivide.
 * Different ChunkSize.
*/
public class MainForKLLBound {
    int dataType=-233;
    static double Pr = 0.80;
    static int startType=0,endType=4;
    static int pageN = 8192;
    static int N = pageN * 6712/*6712*/, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 128; // CHECK IT
    static double[] a;
    static KLLSketchForQuantile[] pageKLL;
    static long[] pageMinV,pageMaxV;
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
                    random.nextDouble();
//                    i-0.5*Math.floor(i/pageN)*pageN;
//                    i;
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
        pageKLL = new KLLSketchForQuantile[pageNum];
        pageMinV = new long[pageNum];pageMaxV = new long[pageNum];

        int enoughMemByte = pageN * 9;
        for (int i = 0; i < pageNum; i++) {
            pageMinV[i] = Long.MAX_VALUE;
            pageMaxV[i] = Long.MIN_VALUE;
            LongKLLSketch worker = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);
            for (int j = 0; j < pageN; j++) {
                long longV =dataToLong(a[i*pageN+j]);
                worker.update(longV);
                pageMinV[i] = Math.min(pageMinV[i],longV);
                pageMaxV[i]=Math.max(pageMaxV[i],longV);
            }
            worker.compactBeforeSerialization();
            pageKLL[i] = worker;
//            worker.show();
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
        double err_full = 0, err_merge_page = 0,bound_merge_page=0;
        double prOutOfBound = 0;
        double[] query_a = new double[queryN];

//        FastKLLSketchLazyForBound.queriedBound.clear();

        for (int T = 0; T < TEST_CASE; T++) {
//            System.out.println("\tTEST_ID:"+T);
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

                int pageL = (L+pageN-1)/pageN, pageR = R/pageN;
                int posL = pageL*pageN, posR = pageR*pageN;
//            System.out.println("\t\t\t"+posL+"\t"+posR);

                merge_page_time -= new Date().getTime();

                FastKLLSketchLazyForBound merge_page_worker = new FastKLLSketchLazyForBound(maxMemoryByte);
//                KLLSketchLazyForBound merge_page_worker = new KLLSketchLazyForBound(maxMemoryByte);

//                for(int i=L;i<Math.min(R,posL);i++)
//                    merge_page_worker.update(dataToLong(a[i]));
//                for (int i=pageL;i<pageR; i++) {
//                    merge_page_worker.mergeWithTempSpace(Collections.singletonList(pageKLL[i]));
//                    merge_page_worker.addRecord(true, pageMinV[i],pageMaxV[i],pageKLL[i].cntLevel-1);
//                }
//                for(int i=Math.max(L,posR);i<R;i++)
//                    merge_page_worker.update(dataToLong(a[i]));


                for(int i=L;i<R;i++)
                    merge_page_worker.update(dataToLong(a[i]));
//                merge_page_worker.show();

                if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
                Arrays.sort(query_a);

                double q_tmp=1e-4,q_start = q_tmp, q_end = 1-q_tmp;
//                q_start=0.2;q_end=0.8;
                double q_add = 1e-4, q_count = Math.floor((q_end - q_start + 1e-10) / q_add) + 1;
                double ratioPerQuery = 1.0/(q_count * TEST_CASE*N/queryN);
                for (double q = q_start; q < q_end + 1e-10; q += q_add) {

                    int query_rank = (int) (q * queryN);

                    double merge_page_v = longToResult(merge_page_worker.findMinValueWithRank(query_rank));
                    int merge_page_delta_rank = getDeltaRank(query_a, queryN, merge_page_v, query_rank);
                    double merge_page_relative_err = 1.0 * merge_page_delta_rank / (queryN);
                    err_merge_page += Math.abs(merge_page_relative_err) * ratioPerQuery;

                    double bound = merge_page_worker.queryBound(dataToLong(merge_page_v),Pr)/queryN;
                    bound_merge_page+=bound* ratioPerQuery;

                    if(Math.abs(merge_page_relative_err)>bound)prOutOfBound+=ratioPerQuery;

//                    System.out.println("Q="+q+"\tapprox v:\t"+merge_page_v+"\t\tbound(ral,rank):\t("+bound+","+bound*queryN+")"+"\t\tactErr(ral,rank):\t("+merge_page_relative_err+","+merge_page_relative_err*queryN+")");

//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
//                    System.out.println("\t\t\trk:"+query_rank+"\t\t"+Math.abs(merge_divide_relative_err));
                }
//                System.out.println("\n");
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("TEST_CASE="+TEST_CASE+"\tDATASET:"+dataType+"\tqueryN=" + queryN + "\tavg_err:" + (int)(err_merge_page*queryN) + "\tavg_bound:"+(int)(bound_merge_page*queryN)+"\t|\tprOutOfBound\t"+prOutOfBound+"\tsetFailPr\t"+(1.0-Pr));
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForKLLBound main;

        System.out.println("TEST for Correct Sketch Divide" + "\n");
        main = new MainForKLLBound();
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            main.RESULT_LINE=0;
            main.prepareA(dataType);
            for (int i : new int[]{/*pageN,chunkN,*/N})
                for (int query_mem : new int[]{/*1024*16,1024*32,1024*64,*/1024*128/*,1024*256*//*,1024*512,1024*1024*/})
                    for (int page_seri : new int[]{/*128,256,512,*/8192/2/*,128,256,512,1024/*,2048,4096,8192*/})
                        main.testError(i, query_mem, page_seri);
                main.RESULT_LINE++;
//            System.out.println("byPage & byChunkDivide\nTEST_CASE=" + TEST_CASE);
        }
//        System.out.println("\nError rate:");
//        for (String s : main.err_result)
//            System.out.println(s);
        System.out.println("\t\t\tALL_TIME:"+(new Date().getTime()-START_T));
    }
}
