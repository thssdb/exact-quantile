import it.unimi.dsi.fastutil.doubles.DoubleArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Random;

/** Compare sketchByPage & sketchByChunkDivide.
 * Different ChunkSize.
*/
public class MainForExactLazyKLLFixedPr {
    int dataType = -233;
    static int startType = 0, endType = 0;
    double fixPr;
    static int pageN = 8192;
    static int N = pageN * 6712, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 2000; // CHECK IT
    static double[] a;
    //    static double[] pageMinV,pageMaxV;
    static int compaction_level;
    ArrayList<String> time_result = new ArrayList<>();
    int RESULT_LINE = 0;
    Random random = new Random(233);

    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] = //longToResult(i&4095);
                    //Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
                    random.nextGaussian();
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

//    public void prepareWorker(int maxSeriesByte) {
//        pageMinV = new double[pageNum];pageMaxV = new double[pageNum];
//
//        int enoughMemByte = pageN * 9;
//        for (int i = 0; i < pageNum; i++) {
//            pageMinV[i] = Long.MAX_VALUE;
//            pageMaxV[i] = Long.MIN_VALUE;
//            for (int j = 0; j < pageN; j++) {
//                double v =a[i*pageN+j];
//                pageMinV[i] = Math.min(pageMinV[i],v);
//                pageMaxV[i]=Math.max(pageMaxV[i],v);
//            }
////            worker.show();
//        }
//    }

    public int getValueActualRank(double[] sortedA, int queryN, double v) { // number of elements <= v
        int L = 0, R = queryN - 1;
        while (L < R) {
            int mid = (L + R + 1) >>> 1;
            if (v < sortedA[mid]) R = mid - 1;
            else L = mid;
        }
        return L;
    }

    public int getValueLessThan(double[] sortedA, int queryN, double v) { // number of elements <= v
        int L = 0, R = queryN - 1;
        while (L < R) {
            int mid = (L + R + 1) >>> 1;
            if (sortedA[mid] < v) L = mid;
            else R = mid - 1;
        }
        return sortedA[L] < v ? L : L - 1;
    }

    public int getDeltaRank(double[] sortedA, int queryN, double v, int targetRank) {
        int rank_L = getValueLessThan(sortedA, queryN, v) + 1;
        int rank_R = getValueActualRank(sortedA, queryN, v);
//        System.out.println("\t\t\t"+targetRank+"\t\tresultLR:"+rank_L+"..."+rank_R+"\t\tresV:"+v);
        if (targetRank >= rank_L && targetRank <= rank_R) return 0;
        else return targetRank < rank_L ? (targetRank - rank_L) : (targetRank - rank_R);
    }

    private long dataToLong(double data) {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }

    private double longToResult(long result) {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }


    public void testError(int queryN, int maxMemoryByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time = 0, merge_page_time = 0;
        double avg_iteration = 0;
//        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        double ALLPrSum = 0, ALLPrCount = 0, FailCount = 0, IterCount = 0;

//        System.out.println("\t!! DataSketchesKLL_K="+DataSketchesKLL_K);
        for (int T = 0; T < TEST_CASE; T++) {
//            System.out.println("\tTEST_ID:"+T);
//            int L = 0, R = queryN;
            int L = LL[T], R = RR[T];

            merge_page_time -= new Date().getTime();
//            KLLSketchLazyExact lazyWorker = new KLLSketchLazyExact(maxMemoryByte);
            KLLSketchLazyExactPriori lazyWorker = new KLLSketchLazyExactPriori(maxMemoryByte);


            for (int i = L; i < R; i++)
                lazyWorker.update(dataToLong(a[i]));
//            lazyWorker.show();

//            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
//            Arrays.sort(query_a);

            double q_tmp = 1e-2, q_start = q_tmp, q_end = 1 - q_tmp;
                q_start=0.50;q_end=0.50;
            double q_add = 1e-2, q_count = Math.floor((q_end - q_start + 1e-10) / q_add) + 1;
            double ratioPerQuery = 1.0 / (q_count * TEST_CASE);
            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
//                if(q<0.189||q>0.194)continue;
//                int query_rank1 = (int) Math.floor(q * queryN),query_rank2 = (int) Math.ceil(q * queryN);
//                int query_rank1 = (int) Math.round(q * queryN),query_rank2 = query_rank1;
                int query_rank1 = (int) Math.floor(q * (queryN - 1) + 1), query_rank2 = (int) Math.ceil(q * (queryN - 1) + 1);
//                System.out.println("\n\t\t\t\t\t\tq rank1,2:"+q+"\t"+query_rank1+","+query_rank2);
                long last_n = queryN;

                double[] deterministic_result, iterate_result = lazyWorker.findResultRange(query_rank1, query_rank2, fixPr);
                deterministic_result = lazyWorker.findResultRange(query_rank1, query_rank2, 1.0);
                avg_iteration += ratioPerQuery;
                IterCount += 1;
                int MMP = 0;
                while (deterministic_result[0] < deterministic_result[1] && deterministic_result.length != 3) {
                    if (++MMP > 50) break;
//                    System.out.println("\n\tnew iter.??\t"+iterate_result[0]+"..."+iterate_result[1]);
                    avg_iteration += ratioPerQuery;
//                    KLLSketchLazyExact cntWorker = new KLLSketchLazyExact(maxMemoryByte);
                    KLLSketchLazyExactPriori cntWorker = new KLLSketchLazyExactPriori(maxMemoryByte);
                    IterCount += 1;
                    double valL = iterate_result[0], valR = iterate_result[1];
                    int CountOfLessThanValL = 0;
                    for (int i = L; i < R; i++) {
                        if (a[i] >= valL && a[i] <= valR) cntWorker.update(dataToLong(a[i]));
                        else if (a[i] < valL) CountOfLessThanValL++;
                    }
//                    cntWorker.show();
//                    System.out.println("\t\t\t\t\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tcntSketch_N:"+cntWorker.getN());
                    int cntRank1 = query_rank1 - CountOfLessThanValL;
                    int cntRank2 = query_rank2 - CountOfLessThanValL;
//                    System.out.println("\t\t\t\t\t\tcntRank:"+cntRank1+" "+cntRank2);
                    if (cntRank1 <= 0 || cntRank2 > cntWorker.getN()) { // iteration failed.
//                        System.out.println("\t\t\t\t\t\titerate fail.");
                        if (cntRank1 <= 0)
                            iterate_result = new double[]{deterministic_result[0], iterate_result[0]};
                        else iterate_result = new double[]{iterate_result[1], deterministic_result[1]};
                        deterministic_result = iterate_result;
                        if (deterministic_result[0] == deterministic_result[1]) continue;
                        FailCount += 1;
                        continue;
                    }
//                    System.out.println("\t\t\t\t\t\titerate success.");
                    iterate_result = cntWorker.findResultRange(cntRank1, cntRank2, fixPr);
                    deterministic_result = cntWorker.findResultRange(cntRank1, cntRank2, 1.0);
                    last_n = cntWorker.getN();
//                    System.out.println("\t\t\t\ttcntL,R:"+iterate_result[0]+","+iterate_result[1]+"\t\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tcntRank:"+cntRank1+" "+cntRank2);
                }
                double exact_quantile_v = (iterate_result[0] + iterate_result[1]) * 0.5;
//                System.out.println("\t\tq:"+q+"\texact_quantile="+exact_quantile_v+"\t\t\titer:"+(MMP+1));
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("fixPR:\t" + fixPr + "\tTEST_CASE=\t" + TEST_CASE + "\tDATASET:\t" + dataType + "\tqueryN:\t" + queryN + "\tmemory:\t" + maxMemoryByte + "\t\tavg_iteration:\t" + avg_iteration + "\t\tavgFailRate:\t" + (FailCount / IterCount));
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForExactLazyKLLFixedPr main;
        DoubleArrayList prList = new DoubleArrayList();
        for (double tmp = 0.50; tmp < 0.7 - 1e-9; tmp += 0.01) {
            prList.add(tmp);
        }
        for (double tmp = 0.70; tmp < 0.99 - 1e-9; tmp += 0.01) {
            prList.add(tmp);
        }
        for (double tmp = 0.99; tmp < 0.999 - 1e-9; tmp += 1e-3) {
            prList.add(tmp);
        }
        for (double tmp = 0.999; tmp < 0.9999 - 1e-9; tmp += 1e-4) {
            prList.add(tmp);
        }
        for (double tmp = 0.9999; tmp <= 0.99999 + 1e-9; tmp += 1e-5) {
            prList.add(tmp);
        }
//        for (double tmp : new double[]{0.5,0.7,0.8,0.85,0.9,0.95,0.99,0.999,0.9999,0.99999}) {
//            prList.add(tmp);
//        }
//        for(double tmp=0.999;tmp<0.9999-1e-9;tmp+=1e-4){
//            prList.add(tmp);
//        }
//        prList.add(1.0);
        main = new MainForExactLazyKLLFixedPr();
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            main.prepareA(dataType);
            for (int queryN : new int[]{100000/*,200000*/})
                for (int query_mem : new int[]{1024*2/*,1024 * 8*/}) {
                    for (double pr : prList) {
                        main.fixPr = pr;
                        main.testError(queryN, query_mem);
                    }
                    System.out.println("\n\n");
                }
            System.out.println("\n-------------------------\n");
        }
        System.out.println("\t\t\tALL_TIME:" + (new Date().getTime() - START_T));
    }
}
