import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;

/** Compare sketchByPage & sketchByChunkDivide.
 * Different ChunkSize.
*/
public class MainForExactLazyKLLFindPrAndCheckWithIID {
    double query_q=0.5;
    int dataType=-233;
    static int startType=1,endType=5;
    static int pageN = 8192;
    static int N = pageN * 6712/*1000000*/, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 32768*4; // CHECK IT
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
                    random.nextDouble();
//                    i-0.5*Math.floor(i/pageN)*pageN;
//                    i;
            return;
        }
        if (dataType == 5) {
            for (int i = 0; i < N; i++)
                a[i] = random.nextGaussian();
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


    public void testError(int queryN, int maxMemoryByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time = 0, merge_page_time = 0;
        double avg_iteration=0;
        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        double ALLPrSum=0,ALLPrCount=0,FailCount=0,IterCount=0;

        int L = 0, R = queryN;

        if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
        Arrays.sort(query_a);
        int query_rank = (int) (query_q * queryN);
        double answer = query_a[query_rank];
//        System.out.println("\t!! DataSketchesKLL_K="+DataSketchesKLL_K);
        IntArrayList rankErrList =new IntArrayList();
        for (int T = 0; T < TEST_CASE; T++) {
//            System.out.println("\tTEST_ID:"+T);
//            int L = 0, R = queryN;

            KLLSketchLazyExact firstWorker = new KLLSketchLazyExact(maxMemoryByte);
            for (int i = L; i < R; i++) firstWorker.update(dataToLong(a[i]));
//            firstWorker.showCompact();

            int answerApproxRank = firstWorker.getApproxRank(dataToLong(answer));
            int answerApproxRankErr = Math.abs(answerApproxRank - query_rank);
//            System.out.println("\t\trankErr:" + answerApproxRankErr + "(" + 1.0 * answerApproxRankErr / queryN + ")");
            rankErrList.add(answerApproxRankErr);
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("CheckErrPR!\tTEST_CASE="+TEST_CASE+"\tDATASET:"+dataType+"\tqueryN:\t" + queryN+"\tmemory:\t" + maxMemoryByte+"\t\tqueryQ:\t"+query_q);
        KLLSketchLazyExact firstWorker = new KLLSketchLazyExact(maxMemoryByte);
        for (int i = L; i < R; i++) firstWorker.update(dataToLong(a[i]));
        rankErrList.sort(Integer::compare);
        DoubleArrayList prList = new DoubleArrayList();
//        for(double pr=0.5;pr<=0.99;pr+=0.05)prList.add(pr);
//        for(double pr=0.96;pr<=0.99+1e-5;pr+=0.01)prList.add(pr);
        for(double pr=0.01;pr<=0.99+1e-5;pr+=0.01)prList.add(pr);
        for(double pr=0.991;pr<=0.999;pr+=0.001)prList.add(pr);
        firstWorker.showCompact();
        for(double pr:prList){
            System.out.println("pr:\t"+pr+"\t\tactualPrErr:\t"+rankErrList.getInt((int)Math.floor(TEST_CASE*(pr)))+"\t\testimated_PrErr:\t"+firstWorker.queryRankErrBound(dataToLong(answer),pr));
        }
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForExactLazyKLLFindPrAndCheckWithIID main;

        main = new MainForExactLazyKLLFindPrAndCheckWithIID();
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            main.RESULT_LINE=0;
            main.prepareA(dataType);
            for (int queryN : new int[]{2000000/*,100000,200000,400000,500000,700000,1000000,2000000,4000000,5000000,7000000,10000000/*,20000000,50000000/**/})
                for (int query_mem : new int[]{1024*8})
                    main.testError(queryN, query_mem);
//                    for(double tmpPr=0.001;tmpPr>=1e-10;tmpPr/=2) {
//                        main.tmpPR=1-tmpPr;
//                        main.testError(queryN, query_mem);
//                    }
                main.RESULT_LINE++;
//            System.out.println("byPage & byChunkDivide\nTEST_CASE=" + TEST_CASE);
        }
//        System.out.println("\nError rate:");
//        for (String s : main.err_result)
//            System.out.println(s);
        System.out.println("\t\t\tALL_TIME:"+(new Date().getTime()-START_T));
    }
}
