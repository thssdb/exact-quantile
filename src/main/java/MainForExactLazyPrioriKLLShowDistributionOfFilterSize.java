import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Date;
import java.util.Random;


public class MainForExactLazyPrioriKLLShowDistributionOfFilterSize {
    int dataType = -233;
    static int startType = 0, endType = 4;
    double fixPr;
    static int pageN = 8192;
    static int N = pageN * 6712, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 64,QUANTILE_PER_TEST=5,Pr_NUM=64; // CHECK IT
    static double errST=0.5,errED=5e-4;
    static double[] a;
    //    static double[] pageMinV,pageMaxV;
    static int compaction_level;
    ObjectArrayList<StringBuilder>RESULT=new ObjectArrayList<>();
    int RESULT_LINE = 0;
    Random random = new Random(233);

    private double getF(long n){
//        return 1.0;
        return n;
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


    public void testError(int queryN, int maxMemoryByte,boolean checkNormal) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        IntArrayList FilterSize=new IntArrayList();
//        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        double ALLPrSum = 0, ALLPrCount = 0, FailCount = 0, prIterCount = 0;

        for (int T = 0; T < TEST_CASE; T++) {
            int L = LL[T], R = RR[T];

//            KLLSketchLazyExact lazyWorker = new KLLSketchLazyExact(maxMemoryByte);
            KLLSketchLazyExactPriori lazyWorker = new KLLSketchLazyExactPriori(maxMemoryByte);


            for (int i = L; i < R; i++)
                lazyWorker.update(dataToLong(a[i]));
//            if(T==0) {
//                lazyWorker.show();
//                System.out.println("\t\t\tmw:\t"+KLLSketchLazyExactPriori.queryRankErrBound(lazyWorker.compactionNumInLevel.toIntArray(),1.0));
//            }

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

                KLLSketchLazyExactPriori cntWorker = new KLLSketchLazyExactPriori(maxMemoryByte);
                if(deterministic_result[0]!=iterate_result[0])prIterCount += 1;
                double valL = iterate_result[0], valR = iterate_result[1];
                int CountOfLessThanValL = 0;
                for (int i = L; i < R; i++) {
                    if (a[i] >/*=*/ valL && a[i] </*=*/ valR) cntWorker.update(dataToLong(a[i]));
                    else if (a[i] < valL) CountOfLessThanValL++;
                }
                FilterSize.add((int)cntWorker.getN());
            }
        }
        double avgFilterSize=0,varFS=0;
        for(int fs:FilterSize)avgFilterSize+=fs;
        avgFilterSize/=FilterSize.size();
        for(int fs:FilterSize)varFS+=Math.pow(fs-avgFilterSize,2.0);
        varFS/=FilterSize.size();

        simuWorker=new KLLSketchLazyEmptyForSimuCompact(maxMemoryByte/8);
        int prERR=KLLSketchLazyExactPriori.queryRankErrBound(simuWorker.simulateCompactNumGivenN(queryN),fixPr);
//        System.out.println("simuLvSize:\t"+ Arrays.toString(simuWorker.lvSize));
        int maxERR=KLLSketchLazyExactPriori.queryRankErrBound(simuWorker.simulateCompactNumGivenN(queryN),1.0);
//        System.out.println("\t\t\tmw:\t"+maxERR);
        double tmpSig2=0;for(int i=0;i<simuWorker.compactNum.length;i++)tmpSig2+=simuWorker.compactNum[i]/2.0*Math.pow(2,i*2);
        int estiFS=prERR*2;

        RESULT.get(RESULT_LINE).append("\t\t"+fixPr+"\t"+avgFilterSize+"\t"+varFS+"\t"+estiFS+"\t"+maxERR+"\t"+tmpSig2+"\t\t\t");
//        System.out.print("\t\t"+fixPr+"\t"+avgFilterSize+"\t"+varFS+"\t"+estiFS+"\t"+maxERR+"\t"+tmpSig2+"\t\t\t");
        if(checkNormal) {
            FilterSize.sort(Integer::compare);
            for (int fs : FilterSize)RESULT.get(RESULT_LINE).append("\t"+fs);
//            for (int fs : FilterSize) System.out.print("\t" + fs);
//            System.out.println();
        }
//        System.out.println();
    }
    KLLSketchLazyEmptyForSimuCompact simuWorker;


    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForExactLazyPrioriKLLShowDistributionOfFilterSize main;
        DoubleArrayList prList = new DoubleArrayList();
        for(double i=0,pNum=Pr_NUM,err=errST,rela=Math.exp(Math.log(errED/errST)/(pNum-1));i<pNum;i++,err*=rela)
            prList.add(1-err);
//        System.out.println("\t\t"+prList);
//        prList.add(1.0);
        main = new MainForExactLazyPrioriKLLShowDistributionOfFilterSize();
        int tmpLine=0,ALL_LINE=prList.size()+8;
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT

            for(int i=0;i<ALL_LINE;i++)main.RESULT.add(new StringBuilder());

            main.prepareA(dataType);
            main.RESULT.get(tmpLine+0).append("\t|||DATASET:"+"\t"+dataType);
            for(int i=1;i<ALL_LINE;i++)main.RESULT.get(tmpLine+i).append("\t|||\t");

            int[] nn=new int[]{100000,400000,1000000};
            int[][] mm=new int[][]{new int[]{1024*7,1024*8,1024*9},new int[]{1024*14,1024*16,1024*18},new int[]{1024*22,1024*25,1024*28}};
            int[] ttcc=new int[]{10,2,1};

            for(int ii=0;ii<3;ii++)
                for(int jj=0;jj<3;jj++){
                    int queryN=nn[ii],queryMem=mm[ii][jj];
                    TEST_CASE*=ttcc[ii];

                    main.RESULT.get(tmpLine+0).append("\t\t\t\t\t\t\t\t\t\t\t");
                    main.RESULT.get(tmpLine+1).append("\t\tqueryN:\t"+queryN+"\tMemory\t"+queryMem+"\t\t\t\t\t");
                    main.RESULT.get(tmpLine+2).append("\t\tPr\tavgFilterSize\tvarFilterSize\testiFS\tsum_mw\tsum_0.5mw^2\t\tdetail\t");
//                    System.out.println("show Var Of FilterSize!\tTEST_CASE=" + TEST_CASE + "\tDATASET:" + dataType + "\tqueryN:\t" + queryN + "\tmemory:\t" + queryMem);
//                    System.out.println("\t\tPr\tavgFilterSize\tvarFilterSize\testiFS\tsum_mw\tsum_0.5mw^2\t\tdetail");

                    main.RESULT_LINE=tmpLine+3;
                    for (double pr : prList) {
                        main.fixPr = pr;
//                        main.testError(queryN, queryMem,false);
                        main.RESULT_LINE++;//break;
                    }
                    if(ii==0&&jj==1)
                    for (double pr : new double[]{0.75,0.85,0.95,0.99}) {
                        main.fixPr = pr;
                        main.RESULT.get(main.RESULT_LINE).append("\tN=\t"+queryN+"\tM=\t"+queryMem+"\t");
                        main.testError(queryN, queryMem,true);
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
