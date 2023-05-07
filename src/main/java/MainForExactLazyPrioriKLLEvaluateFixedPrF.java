import it.unimi.dsi.fastutil.doubles.DoubleArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.Random;


public class MainForExactLazyPrioriKLLEvaluateFixedPrF {
    static boolean DEBUG_PRINT = false;
    int dataType = -233;
    static int startType = 0, endType = 0;
    static int pageN = 8192;
    static int N = /*81000000*/8192 * 6712, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 1,Pr_NUM=64; // CHECK IT
    static double errST=0.5,errED=5e-4;
    static double[] a;
    static DoubleArrayList prList = null;
    ArrayList<String> time_result = new ArrayList<>();
    int RESULT_LINE = 0;
    Random random = new Random(233);

    static int DEBUG_COUNT_SIMULATION = 0;

    private double getF(long n){
        return 1.0;
//        return n;
    }

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
        if (dataType == 4) {
            for (int i = 0; i < N; i++)
                a[i] = i;
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
//
//    KLLSketchLazyEmptyForSimuCompact simuWorker;
//
////    private double[] simulateIteration(double casePr, double lastPr, int depth, int n, int maxMemoryNum, int[] compactNum) {
////        if (n <= 0) return new double[]{0, 1.0};
////        if (n <= maxMemoryNum) return
//////            new double[]{
////////            Math.max(0.75, Math.log(n) / Math.log(maxMemoryNum)),
//////                1.0,
//////                1.0};
////            new double[]{
//////            Math.max(0.75, Math.log(n) / Math.log(maxMemoryNum)),
////                getF(n),1.0
////            };
////        int maxERR = 0;
////        for (int i = 0; i < compactNum.length; i++) maxERR += compactNum[i] << i;
////        double bestSimuIter = 1e3, bestSimuPr = 1.0;
////
////        double pr = bestSimuPr = lastPr;
////        int prERR = KLLSketchLazyExactPriori.queryRankErrBound(compactNum, pr);
////        int successN = prERR * 2;
////        int failN = (Math.min(n, maxERR * 2) - prERR) / 2;
//////            KLLSketchLazyEmptyForSimuCompact simuWorker = new KLLSketchLazyEmptyForSimuCompact(/*n, */maxMemoryNum);
////        int[] succComNum, failComNum;
////        if (failN < successN) {
////            failComNum = simuWorker.simulateCompactNumGivenN(failN);
////            succComNum = simuWorker.simulateCompactNumGivenN(successN);
////        } else {
////            succComNum = simuWorker.simulateCompactNumGivenN(successN);
////            failComNum = simuWorker.simulateCompactNumGivenN(failN);
////        }
//////        bestSimuIter = 1 + pr * simulateIteration(casePr * pr, pr, depth + 1, successN, maxMemoryNum, succComNum)[0] + (1 - pr) * (1 + simulateIteration(casePr * (1 - pr), pr, depth + 1, failN, maxMemoryNum, failComNum)[0]);
////        bestSimuIter = getF(n) + pr * simulateIteration(casePr * pr, pr, depth + 1, successN, maxMemoryNum, succComNum)[0] + (1 - pr) * (getF(successN) + simulateIteration(casePr * (1 - pr), pr, depth + 1, failN, maxMemoryNum, failComNum)[0]);
////
////        DEBUG_COUNT_SIMULATION++;
////        return new double[]{bestSimuIter, bestSimuPr};
////    }
//    private double simulateIteration(double casePr, double fixPr, int depth, int maxMemoryNum,int avgN,double sig2,long maxERR) {
//        if (avgN <= 0) return 0;
//        if (avgN + maxERR <= maxMemoryNum) return 1.0;
//        NormalDistribution normalDis = NormalDistribution.of(avgN, Math.sqrt(sig2));
//        double finishPr = normalDis.cumulativeProbability(maxMemoryNum);
////        System.out.println("\t\t\t\t?????\t\tavgN:"+avgN+"\t\tmaxMem:"+maxMemoryNum+"\t\tfinishPr:"+finishPr);
//        double continuePr = 1 - finishPr;
//        if (continuePr < 1e-5) return 1.0;
//        TruncatedNormalDistribution tnd = TruncatedNormalDistribution.of(avgN, Math.sqrt(sig2), maxMemoryNum, Double.MAX_VALUE);
//        int continueAvgN = (int) Math.ceil(tnd.getMean());
//
//        int[] conCompactNum = simuWorker.simulateCompactNumGivenN(continueAvgN);
//        double conSig2 = KLLSketchLazyEmptyForSimuCompact.getSig2(conCompactNum);
//        long conMaxErr = KLLSketchLazyEmptyForSimuCompact.getMaxError(conCompactNum);
//
////        System.out.println("\tdepth:"+depth+" avgN:"+avgN+"\tfixPr:"+fixPr+"\t\tcontinuePr:"+continuePr+"\t\t\tconAvgN:"+continueAvgN+"\t\tconMaxERR:"+conMaxErr);
//
//
//        int conSuccessN = conMaxErr==0?1:KLLSketchLazyExactPriori.queryRankErrBoundGivenParameter(conSig2, conMaxErr, fixPr) * 2;
//        int conFailN = (Math.min(continueAvgN, (int) conMaxErr * 2) - conSuccessN) / 2;
////        System.out.println("\t\t\tconSuccessN:"+conSuccessN+"\t\tconFailN:"+conFailN);
////        bestSimuIter = 1 + pr * simulateIteration(casePr * pr, pr, depth + 1, successN, maxMemoryNum, succComNum)[0] + (1 - pr) * (1 + simulateIteration(casePr * (1 - pr), pr, depth + 1, failN, maxMemoryNum, failComNum)[0]);
//        double simuPrF =
//            finishPr * 1.0 +
//                continuePr * (
//                    getF(continueAvgN) + fixPr * simulateIteration(casePr * fixPr, fixPr, depth + 1, maxMemoryNum, conSuccessN, conSig2, conMaxErr)
//                        + (1 - fixPr) * (getF(conSuccessN) + simulateIteration(casePr * (1 - fixPr), fixPr, depth + 1, maxMemoryNum, conFailN, conSig2, conMaxErr))
//                );
//        DEBUG_COUNT_SIMULATION++;
//        return simuPrF;
//    }
//
//    private double evaluatePr(int maxMemoryNum, double fixPr, int detN) {
//        if(detN<=maxMemoryNum)return 1.0;
//        int maxERR = KLLSketchLazyExactPriori.queryRankErrBound(simuWorker.simulateCompactNumGivenN(detN), 1.0);
//        double sig2 = simuWorker.getSig2();
//        long nextMaxErr=simuWorker.getMaxError();
//        int prERR = KLLSketchLazyExactPriori.queryRankErrBound(simuWorker.compactNum, fixPr);
//
//        int successAvgN = prERR * 2;
//        int failAvgN = (Math.min(detN, maxERR * 2) - successAvgN) / 2;
//
////        System.out.println("\t\tfirst paas. successAvgN:\t"+successAvgN+"\t\tfailAvgN:\t"+failAvgN);
////        System.out.println("\t\t\t\t\t"+sig2+"\t"+maxERR+"\t\t\tnextMaxERR:\t"+nextMaxErr);
//
//        double simulateResult;
//        double successResult = simulateIteration(fixPr, fixPr, 1, maxMemoryNum, successAvgN,sig2,nextMaxErr);
//        double failResult = simulateIteration((1 - fixPr), fixPr, 1, maxMemoryNum, failAvgN, sig2, nextMaxErr);
//        simulateResult = getF(detN) + fixPr * successResult + (1 - fixPr) * (getF(successAvgN) + failResult);
////        if(DEBUG_PRINT)System.out.println("\t\t\t\t\t\t\t\tcntPR:"+fixPr+"\tsuccessN:\t"+succN+"\t\tfailN:\t"+failN+/*"\t\testi_iter:\t"+estimateIterationNum+*/"\t\tsimu_iter:\t"+simulateResult[0]+"\tsimu_nextSuccessPr:"+simulateResult[1]);
//        return simulateResult;
//    }



    public void testError(int queryN, int maxMemoryByte,int merge_buffer_ratio,int multi_quantile) throws IOException {
//        System.out.println("FindFixPR!\tTEST_CASE=" + TEST_CASE + "\tDATASET:" + dataType + "\tqueryN:\t" + queryN + "\tmemory:\t" + maxMemoryByte);
//        System.out.println("\t\tPr\testiF(N,Mem,Pr)");
////
//        double avgF = 0,avgPass=0;
//        double[] query_a = new double[queryN];
//        simuWorker=new KLLSketchLazyEmptyForSimuCompact(maxMemoryByte/8);
//        for(double fixPr:prList) {
//            double estiPass = PrioriBestPrHelper.evaluatePrFromScratch(maxMemoryByte/8,fixPr,queryN,merge_buffer_ratio,multi_quantile);
//            System.out.println("\t\t" + fixPr+ "\t" + estiPass );
//        }
        PrioriBestPrHelper helper = new PrioriBestPrHelper(maxMemoryByte/8,queryN,merge_buffer_ratio,multi_quantile);
        double[] findPrResult = helper.findBestPr(5e-4,5e-4,5e-1,queryN);
        System.out.println("MULTI_QUANTILE:\t"+multi_quantile+"\tbestPr:\t" + findPrResult[0]+"\t"+findPrResult[1] );
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        prList = new DoubleArrayList();

        for(double i=0,pNum=Pr_NUM,err=errST,rela=Math.exp(Math.log(errED/errST)/(pNum-1));i<pNum;i++,err*=rela) {
            prList.add(1 - err);
//            System.out.println("\t\t\t+\t+delta:\t"+err+"\t\t"+Math.log(err/errED)*(pNum-1)/Math.log(errST/errED)+"\t\t"+Math.log(errST/errED));
        }

        long START_T = new Date().getTime();
        MainForExactLazyPrioriKLLEvaluateFixedPrF main;

        main = new MainForExactLazyPrioriKLLEvaluateFixedPrF();

        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
//            main.prepareA(dataType);
            // 1e7 100quantiles 265KB           1e7 1600Quantiles 4MB
            // 1e6->5e6 1quantile  32KB  best_passes: 2->3
            // multi-peak: 1e6 14KB or 15KB
            int[] nn=new int[]{100000,400000,1000000000};
            int[][] mm=new int[][]{new int[]{1024*7,1024*8,1024*9},new int[]{1024*14,1024*16,1024*18},new int[]{1024*1024*30}};

            for(int ii=2;ii<3;ii++)
                for(int jj=0;jj<1;jj++){
                    int queryN=nn[ii],queryMem=mm[ii][jj];

                    System.out.println("\t\tqueryN:\t"+queryN+"\tMemory\t"+queryMem+"\t");
                    System.out.println("\t\tPr\testiPass");
//                    for(int i=0;i<5;i++)
                    main.testError(queryN, queryMem,5,800);
                    main.testError(queryN, queryMem,0,800);
//                    System.out.println("\t\tDEBUG_COUNT_SIMULATION:\t"+DEBUG_COUNT_SIMULATION);
                }

//            System.out.println("\n-------------------------\n");
        }
        System.out.println("\t\t\tALL_TIME:" + (new Date().getTime() - START_T));
    }
}
