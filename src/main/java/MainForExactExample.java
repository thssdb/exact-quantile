import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;


public class MainForExactExample {
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
        a=new double[N];
        for (int i = 0; i < N; i++)
            a[i] = i;
        XoRoShiRo128PlusRandom random = new XoRoShiRo128PlusRandom();
        for(int i=1;i<N;i++) {
            int j = random.nextInt(i);
            double tmp = a[j];
            a[j] = a[i];
            a[i] = tmp;
        }
//        System.out.println("A:\t"+ Arrays.toString(a));
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
    public void testSample(int queryN, int maxMemoryByte,int merge_buffer_ratio,int multi_quantile){
        PrioriBestPrHelper helper = new PrioriBestPrHelper(maxMemoryByte/8,queryN,merge_buffer_ratio,multi_quantile);
        double[] findPrResult = helper.findBestPr(5e-4,5e-4,5e-1,queryN);
        System.out.println("MULTI_QUANTILE:\t"+multi_quantile+"\tbestPr:\t" + findPrResult[0]+"\t"+findPrResult[1] );

        KLLSketchLazyExactPriori sketch = new KLLSketchLazyExactPriori(maxMemoryByte);
        for(int i=0;i<queryN;i++)sketch.update(dataToLong(a[i]));
        sketch.show();
        sketch.showCompact();

        System.out.println("\t\t\tsketch VarErr:\t"+ KLLSketchLazyEmptyForSimuCompact.getSig2(sketch.getRelatedCompactNum()));
        double q=0.5;
        int query_rank1 = (int) Math.floor(q * (queryN - 1) + 1), query_rank2 = (int) Math.ceil(q * (queryN - 1) + 1);

        double bestPr=findPrResult[0];
        double[] pr_result = sketch.getFilter(0,0,0,0,query_rank1,query_rank2,bestPr);
        double[] pr_result_99 = sketch.getFilter(0,0,0,0,query_rank1,query_rank2,0.99);
        double[] deterministic_result = sketch.getFilter(0,0,0,0,query_rank1,query_rank2,1.0);

        System.out.println("\t\tpr_result:\t"+ Arrays.toString(pr_result));
        System.out.println("\t\t0.99_result:\t"+ Arrays.toString(pr_result_99));
        System.out.println("\t\tdet_result:\t"+ Arrays.toString(deterministic_result));

//        for(double fixPr:prList) {
//            double estiPass = PrioriBestPrHelper.evaluatePrFromScratch(maxMemoryByte/8,fixPr,queryN,merge_buffer_ratio,multi_quantile);
//            System.out.println("\t\t" + fixPr+ "\t" + estiPass );
//        }
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
        MainForExactExample main;

        main = new MainForExactExample();

        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
//            main.prepareA(dataType);
            // 1e7 100quantiles 265KB           1e7 1600Quantiles 4MB
            // 1e6->5e6 1quantile  32KB  best_passes: 2->3
            // multi-peak: 1e6 14KB or 15KB
            int[] nn=new int[]{0,0,200000};
            int[][] mm=new int[][]{new int[]{1024*7,1024*8,1024*9},new int[]{1024*14,1024*16,1024*18},new int[]{1024*8}};

            for(int ii=2;ii<3;ii++)
                for(int jj=0;jj<1;jj++){
                    int queryN=nn[ii],queryMem=mm[ii][jj];
                    N=queryN;
                    main.prepareA(0);
                    System.out.println("\t\tqueryN:\t"+queryN+"\tMemory\t"+queryMem+"Byte\t");
//                    for(int i=0;i<5;i++)
                    main.testSample(queryN, queryMem,0,1);
//                    System.out.println("\t\tDEBUG_COUNT_SIMULATION:\t"+DEBUG_COUNT_SIMULATION);
                }

//            System.out.println("\n-------------------------\n");
        }
        System.out.println("\t\t\tALL_TIME:" + (new Date().getTime() - START_T));
    }
}
