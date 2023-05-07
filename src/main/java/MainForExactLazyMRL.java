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
public class MainForExactLazyMRL {
    int dataType=-233;
    static int startType=0,endType=0;
    static int pageN = 8192;
    static int N = pageN * 6712/*6712*/, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 10,QUANTILE_PER_TEST=5; // CHECK IT
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
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte);
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

        for (int T = 0; T < TEST_CASE; T++) {
//            System.out.println("\tTEST_ID:"+T);
//            int L = 0, R = queryN;
            int L=LL[T],R=RR[T];

            merge_page_time -= new Date().getTime();
            MRLSketchLazy lazyWorker = new MRLSketchLazy(maxMemoryByte);

            for (int i = L; i < R; i++)
                lazyWorker.update(dataToLong(a[i]));
//            lazyWorker.show();
//            lazyWorker.showCompact();

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);

            int q_count=QUANTILE_PER_TEST;
            double ratioPerQuery = 1.0 / (q_count * TEST_CASE);
            for(int INNER_T=0;INNER_T<q_count;INNER_T++)
            {
                double q = random.nextDouble();
                int query_rank1 = (int) Math.floor(q * (queryN-1)+1),query_rank2 = (int) Math.ceil(q * (queryN-1)+1);
//                System.out.println("\n\t\t\t\t\t\tq rank1,2:"+q+"\t"+query_rank1+","+query_rank2);

//                int MMM0=0,MMM1=0,MMMx=0;
//                for(int i=L;i<R;i++){
//                    if(a[i]<=29.87877655)MMM0++;
//                    if(a[i]<=30)MMM1++;
//                    if(a[i]<30)MMMx++;
//                }System.out.println("\t\tMMM0,1:"+MMM0+"\t"+MMM1+"\t\t"+MMMx);
                long last_n=queryN;

                double[] deterministic_result;
                deterministic_result = lazyWorker.findResultRange(query_rank1,query_rank2);
                avg_iteration+=ratioPerQuery;
                int MMP=0;
                while(deterministic_result[0]<deterministic_result[1]&&deterministic_result.length!=3) {
                    if(++MMP>50){
                        System.out.println("\t\t\t\t\t\titerate fail.\titer:"+MMP+"\t\tq:"+q);
                        return;
                    }
//                    System.out.println("\n\tnew iter.??\t"+deterministic_result[0]+"..."+deterministic_result[1]+"\t\tFilterSize:\t"+last_n);
                    avg_iteration+=ratioPerQuery;
                    MRLSketchLazy cntWorker = new MRLSketchLazy(maxMemoryByte);
                    double valL=deterministic_result[0],valR=deterministic_result[1];
                    int CountOfLessThanValL = 0,CountOfValL=0,CountOfValR=0;
                    for(int i=L;i<R;i++) {
                        if (a[i] > valL && a[i] < valR) cntWorker.update(dataToLong(a[i]));
                        else if (a[i] < valL) CountOfLessThanValL++;
                        else if(a[i]==valL)CountOfValL++;
                        else if(a[i]==valR)CountOfValR++;
                    }
//                    cntWorker.show();
//                    System.out.println("\t\t\t\t\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tcntSketch_N:"+cntWorker.getN());
                    int cntRank1=query_rank1-CountOfLessThanValL;
                    int cntRank2=query_rank2-CountOfLessThanValL;
//                    System.out.println("\t\t\t\t\t\tcntRank:"+cntRank1+" "+cntRank2);
                    if(cntRank1<=0||cntRank2>CountOfValL+CountOfValR+cntWorker.getN()){ // iteration failed.
                        System.out.println("\t\t\t\t\t\titerate fail.");
                        return;
                    }
//                    System.out.println("\t\t\t\t\t\titerate success.");
//                    deterministic_result = cntWorker.findResultRange(cntRank1,cntRank2);
                    deterministic_result = cntWorker.getFilter(CountOfValL,CountOfValR,valL,valR,cntRank1,cntRank2);
                    last_n=cntWorker.getN();
//                    System.out.println("\t\t\t\ttcntL,R:"+deterministic_result[0]+","+deterministic_result[1]+"\t\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tcntRank:"+cntRank1+" "+cntRank2+"\t\tcntIterN:"+cntWorker.getN()+"\t\tcntIterErr:"+cntWorker.queryErrBound());
                }
                double exact_quantile_v = (deterministic_result[0]+deterministic_result[1])*0.5;
                System.out.println("\t\tq:"+q+"\texact_quantile="+exact_quantile_v+"\t\t\titer:"+(MMP+1));
            }
        }
        System.out.println("\t\t"+avg_iteration);
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForExactLazyMRL main;

        main = new MainForExactLazyMRL();
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            System.out.println("\t|||DATASET:"+"\t"+dataType);
            main.RESULT_LINE=0;
            main.prepareA(dataType);
            int[] nn=new int[]{100000,400000,1000000};
            int[][] mm=new int[][]{new int[]{1024*7,1024*8,1024*9},new int[]{1024*14,1024*16,1024*18},new int[]{1024*64}};
            int[] ttcc=new int[]{10,2,1};

            for(int ii=2;ii<3;ii++)
                for(int jj=0;jj<1;jj++){
                    int queryN=nn[ii],queryMem=mm[ii][jj];
                    TEST_CASE*=ttcc[ii];
                    System.out.println("\t\tqueryN:\t"+queryN+"\tMemory\t"+queryMem+"\t");
                    System.out.println("\t\tavgPass");

                    main.testError(queryN, queryMem);
                    TEST_CASE/=ttcc[ii];
                }
        }
        System.out.println("\t\t\tALL_TIME:"+(new Date().getTime()-START_T));
    }
}
