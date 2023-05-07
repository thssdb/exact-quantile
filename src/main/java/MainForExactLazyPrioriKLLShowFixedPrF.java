import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Date;
import java.util.Random;


public class MainForExactLazyPrioriKLLShowFixedPrF {
    int dataType = -233;
    static int startType = 0, endType = 3;
    double fixPr;
    static int pageN = 8192;
    static int N = pageN * 6712, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 16*40*5,QUANTILE_PER_TEST=5,Pr_NUM=64; // CHECK IT
    static double errST=0.5,errED=5e-4;
    static double[] a;
    //    static double[] pageMinV,pageMaxV;
    static int compaction_level;
    ObjectArrayList<StringBuilder>RESULT=new ObjectArrayList<>();
    int RESULT_LINE = 0;
    Random random = new Random(233);

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


    public void testError(int queryN, int maxMemoryByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time = 0, merge_page_time = 0;
        double avgF = 0,avgPass=0;
//        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        double ALLPrSum = 0, ALLPrCount = 0, FailCount = 0, prIterCount = 0;

//        System.out.println("\t!! DataSketchesKLL_K="+DataSketchesKLL_K);
        for (int T = 0; T < TEST_CASE; T++) {
            int L = LL[T], R = RR[T];

            merge_page_time -= new Date().getTime();
//            KLLSketchLazyExact lazyWorker = new KLLSketchLazyExact(maxMemoryByte);
            KLLSketchLazyExactPriori lazyWorker = new KLLSketchLazyExactPriori(maxMemoryByte);


            for (int i = L; i < R; i++)
                lazyWorker.update(dataToLong(a[i]));
//            lazyWorker.show();

            int q_count=QUANTILE_PER_TEST;
            double ratioPerQuery = 1.0 / (q_count * TEST_CASE);
            for(int INNER_T=0;INNER_T<q_count;INNER_T++)
            {
                double q = random.nextDouble();

                int query_rank1 = (int) Math.floor(q * (queryN - 1) + 1), query_rank2 = (int) Math.ceil(q * (queryN - 1) + 1);
//                System.out.println("\n\t\t\t\t\t\tq rank1,2:"+q+"\t"+query_rank1+","+query_rank2);
                long last_n = queryN;

                double[] deterministic_result, iterate_result = lazyWorker.findResultRange(query_rank1, query_rank2, fixPr);
                deterministic_result = lazyWorker.findResultRange(query_rank1, query_rank2, 1.0);
                avgF += getF(queryN)*ratioPerQuery;
                avgPass+=ratioPerQuery;
//                prIterCount += 1;
                int MMP = 0;
                while (deterministic_result[0] < deterministic_result[1] && deterministic_result.length != 3) {
                    if (++MMP > 50) break;
//                    System.out.println("\n\tnew iter.??\t"+iterate_result[0]+"..."+iterate_result[1]);
//                    KLLSketchLazyExact cntWorker = new KLLSketchLazyExact(maxMemoryByte);
                    KLLSketchLazyExactPriori cntWorker = new KLLSketchLazyExactPriori(maxMemoryByte);
                    if(deterministic_result[0]!=iterate_result[0])prIterCount += 1;
                    double valL = iterate_result[0], valR = iterate_result[1];
                    int CountOfLessThanValL = 0,CountOfValL=0,CountOfValR=0;
                    for (int i = L; i < R; i++) {
                        if (a[i] > valL && a[i] < valR) cntWorker.update(dataToLong(a[i]));
                        else if (a[i] < valL) CountOfLessThanValL++;
                        else if(a[i]==valL)CountOfValL++;
                        else if(a[i]==valR)CountOfValR++;
                    }
                    avgF += getF(cntWorker.getN())*ratioPerQuery;
                    avgPass+=ratioPerQuery;
//                    cntWorker.show();
//                    System.out.println("\t\t\t\t\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tcntSketch_N:"+cntWorker.getN());
                    int cntRank1 = query_rank1 - CountOfLessThanValL;
                    int cntRank2 = query_rank2 - CountOfLessThanValL;
//                    System.out.println("\t\t\t\t\t\tcntRank:"+cntRank1+" "+cntRank2);
                    if (cntRank1 <= 0 || cntRank2 > CountOfValL+CountOfValR+cntWorker.getN()) { // iteration failed.
//                        System.out.println("\t\t\t\t\t\titerate fail.");
                        FailCount += 1;
                        if (cntRank1 <= 0)
                            iterate_result = new double[]{deterministic_result[0], iterate_result[0]};
                        else iterate_result = new double[]{iterate_result[1], deterministic_result[1]};
                        deterministic_result = iterate_result;
                        if (deterministic_result[0] == deterministic_result[1]) continue;
                        continue;
                    }
//                    System.out.println("\t\t\t\t\t\titerate success.");
                    iterate_result = cntWorker.getFilter(CountOfValL,CountOfValR,valL,valR,cntRank1,cntRank2,fixPr);
                    deterministic_result = cntWorker.getFilter(CountOfValL,CountOfValR,valL,valR,cntRank1,cntRank2,1.0);
//                    System.out.println("\t\t\t\ttcntL,R:"+iterate_result[0]+","+iterate_result[1]+"\tdetL,R:"+deterministic_result[0]+","+deterministic_result[1]+"\t\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tcntRank:"+cntRank1+" "+cntRank2);

//                    iterate_result = cntWorker.findResultRange(cntRank1, cntRank2, fixPr);
//                    deterministic_result = cntWorker.findResultRange(cntRank1, cntRank2, 1.0);
//                    System.out.println("\t\t\t\ttcntL,R:"+iterate_result[0]+","+iterate_result[1]+"\t\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tcntRank:"+cntRank1+" "+cntRank2);
                }
                double exact_quantile_v = (iterate_result[0] + iterate_result[1]) * 0.5;
//                System.out.println("\t\tq:"+q+"\texact_quantile="+exact_quantile_v+"\t\t\titer:"+(MMP+1));
            }
        }
        simuWorker=new KLLSketchLazyEmptyForSimuCompact(maxMemoryByte/8);
        simuWorker.simulateCompactNumGivenN(queryN);
//        System.out.println(Arrays.toString(simuWorker.compactNum));
        int maxERR=KLLSketchLazyExactPriori.queryRankErrBound(simuWorker.compactNum,1.0);
        int prERR=KLLSketchLazyExactPriori.queryRankErrBound(simuWorker.compactNum,fixPr);
        int successN = prERR * 2;
        int failN = (Math.min(queryN, maxERR * 2) - prERR) / 2;

        int[] sucCom= simuWorker.simulateCompactNumGivenN(successN);
        int[] failCom= simuWorker.simulateCompactNumGivenN(failN);
        double estiPrF=getF(queryN)+evaluatePr(maxMemoryByte/8,fixPr,successN,failN,sucCom,failCom)[0];
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        RESULT.get(RESULT_LINE).append("\t\t"+avgF+"\t"+avgPass+"\t"+estiPrF+"\t"+(FailCount / prIterCount)+"\t"+(1-fixPr));
        System.out.println("\t\t"+avgF+"\t"+avgPass+"\t"+estiPrF+"\t"+(FailCount / prIterCount)+"\t"+(1-fixPr));
//        System.out.println("fixPR:\t" + fixPr + "\tTEST_CASE=\t" + TEST_CASE + "\tDATASET:\t" + dataType + "\tqueryN:\t" + queryN + "\tmemory:\t" + maxMemoryByte + "\t\tavgF:\t" + avgF + "\t\tavgFailRate:\t" + (FailCount / prIterCount));
    }

    KLLSketchLazyEmptyForSimuCompact simuWorker;

    private double[] simulateIteration(double casePr, double lastPr, int depth, int n, int maxMemoryNum, int[] compactNum) {
        if (n <= 0) return new double[]{0, 1.0};
        if (n <= maxMemoryNum) return
//            new double[]{
////            Math.max(0.75, Math.log(n) / Math.log(maxMemoryNum)),
//                1.0,
//                1.0};
            new double[]{
//            Math.max(0.75, Math.log(n) / Math.log(maxMemoryNum)),
                getF(n),1.0
            };
        int maxERR = 0;
        for (int i = 0; i < compactNum.length; i++) maxERR += compactNum[i] << i;
        double bestSimuIter = 1e3, bestSimuPr = 1.0;

        double pr = bestSimuPr = lastPr;
        int prERR = KLLSketchLazyExactPriori.queryRankErrBound(compactNum, pr);
        int successN = prERR * 2;
        int failN = (Math.min(n, maxERR * 2) - prERR) / 2;
//            KLLSketchLazyEmptyForSimuCompact simuWorker = new KLLSketchLazyEmptyForSimuCompact(/*n, */maxMemoryNum);
        int[] succComNum, failComNum;
        if (failN < successN) {
            failComNum = simuWorker.simulateCompactNumGivenN(failN);
            succComNum = simuWorker.simulateCompactNumGivenN(successN);
        } else {
            succComNum = simuWorker.simulateCompactNumGivenN(successN);
            failComNum = simuWorker.simulateCompactNumGivenN(failN);
        }
//        bestSimuIter = 1 + pr * simulateIteration(casePr * pr, pr, depth + 1, successN, maxMemoryNum, succComNum)[0] + (1 - pr) * (1 + simulateIteration(casePr * (1 - pr), pr, depth + 1, failN, maxMemoryNum, failComNum)[0]);
        bestSimuIter = getF(n) + pr * simulateIteration(casePr * pr, pr, depth + 1, successN, maxMemoryNum, succComNum)[0] + (1 - pr) * (getF(successN) + simulateIteration(casePr * (1 - pr), pr, depth + 1, failN, maxMemoryNum, failComNum)[0]);

        return new double[]{bestSimuIter, bestSimuPr};
    }


    private double[] evaluatePr(int maxMemoryNum, double Pr, int succN, int failN, int[] succComNum, int[] failComNum) {
        double[] simulateResult = new double[3];
        double[] successResult = simulateIteration(Pr, Pr, 1, succN, maxMemoryNum, succComNum);
        double[] failResult = simulateIteration(0 * (1 - Pr), Pr, 1, failN, maxMemoryNum, failComNum);
        simulateResult[0] = Pr * successResult[0] + (1 - Pr) * (getF(succN) + failResult[0]);
        simulateResult[1] = successResult[1];
        simulateResult[2] = failResult[1];
//        if(DEBUG_PRINT)System.out.println("\t\t\t\t\t\t\t\tcntPR:"+Pr+"\tsuccessN:\t"+succN+"\t\tfailN:\t"+failN+/*"\t\testi_iter:\t"+estimateIterationNum+*/"\t\tsimu_iter:\t"+simulateResult[0]+"\tsimu_nextSuccessPr:"+simulateResult[1]);
        return simulateResult;
    }



    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForExactLazyPrioriKLLShowFixedPrF main;
        DoubleArrayList prList = new DoubleArrayList();
        for(double i=0,pNum=Pr_NUM,err=errST,rela=Math.exp(Math.log(errED/errST)/(pNum-1));i<pNum;i++,err*=rela)
            prList.add(1-err);

//        prList.add(1.0);
        main = new MainForExactLazyPrioriKLLShowFixedPrF();
        int tmpLine=0,ALL_LINE=prList.size()+4;
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT

            for(int i=0;i<ALL_LINE;i++)main.RESULT.add(new StringBuilder());

            main.prepareA(dataType);
            main.RESULT.get(tmpLine+0).append("\t|||DATASET:"+"\t"+dataType);
            for(int i=1;i<ALL_LINE;i++)main.RESULT.get(tmpLine+i).append("\t|||\t");


            int[] nn=new int[]{100000,400000,1000000,1000000};
            int[][] mm=new int[][]{new int[]{1024*8,1024*8,1024*9},new int[]{1024*16,1024*16,1024*18},new int[]{1024*22,1024*28},new int[]{1024*28}};
            int[] ttcc=new int[]{10,3,1,1};

            for(int ii=0;ii<1;ii++)
                for(int jj=0;jj<1;jj++){
                    int queryN=nn[ii],queryMem=mm[ii][jj];
                    TEST_CASE*=ttcc[ii];

                    main.RESULT.get(tmpLine+0).append("\t\t\t\t\t\t\t\t\t\t\t");
                    main.RESULT.get(tmpLine+1).append("\t\tqueryN:\t"+queryN+"\tMemory:\t"+queryMem/1024+"KB\t");
                    main.RESULT.get(tmpLine+2).append("\t\tactualPrF\tactualPass\testiPrF\tavgFailRate\tdelta");
                    System.out.println("\t\tqueryN:\t"+queryN+"\tMemory\t"+queryMem+"\t");
                    System.out.println("\t\tactualPrF\tactualPass\testiPrF\tavgFailRate\tdelta");
//                    System.out.println("show Var Of FilterSize!\tTEST_CASE=" + TEST_CASE + "\tDATASET:" + dataType + "\tqueryN:\t" + queryN + "\tmemory:\t" + queryMem);
//                    System.out.println("\t\tPr\tavgFilterSize\tvarFilterSize\testiFS\tsum_mw\tsum_0.5mw^2\t\tdetail");

                    main.RESULT_LINE=tmpLine+3;
                    for(int i=0;i<Pr_NUM;i++){
                        if(i>=18)break;
                        main.fixPr = prList.getDouble(Pr_NUM-i-1);
                        main.testError(queryN, queryMem);
                        main.RESULT_LINE++;//break;
                    }
                    TEST_CASE/=ttcc[ii];
                }

            tmpLine+=ALL_LINE;
//            System.out.println("\n-------------------------\n");
        }

        for(StringBuilder sb:main.RESULT)System.out.println("\t"+sb.toString());
        System.out.println("\t\t\tALL_TIME:" + (new Date().getTime() - START_T));
    }
}