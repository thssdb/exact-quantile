import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;

/** Compare sketchByPage & sketchByChunkDivide.
 * Different ChunkSize.
*/
public class MainForExactQuickSelectMultiQuantiles {
    static int MULTI_QUANTILES;
    static boolean DEBUG_PRINT=false;
    int dataType=-233;
    static int startType=5,endType=5;
    static int pageN = 4096,batchN=4096;
    static int N = /*81000000*/8192*6713*3, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 1; // CHECK IT
    static double[] a,query_q,answer_v;
    static int[] query_rank1,query_rank2;
    LongArrayList data;
    ArrayList<String> time_result = new ArrayList<>();
    int RESULT_LINE = 0;
    public XoRoShiRo128PlusRandom random = new XoRoShiRo128PlusRandom(233);

    static int DEBUG_COUNT_SIMULATION=0;

    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] = //longToResult(i&4095);
                    //Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
//                    random.nextGaussian();
//                    i-0.5*Math.floor(i/pageN)*pageN;
                    i;
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

    private long dataToLong(double data) {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }
    private double longToResult(long result) {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }


    public void getKth(int L, int R, int qL,int qR) {
        if(L==R){
            long v = data.getLong(L);
            for(int i=qL;i<qR;i++){
                if(query_rank1[i]>=L)answer_v[i]+=longToResult(v);
                if(query_rank2[i]<R)answer_v[i]+=longToResult(v);
            }
            return;
        }
        int pos = L + random.nextInt(R - L);
        long pivot_v = data.getLong(pos), swap_v;

        int leP = L, eqR = R;
        data.set(pos, data.set(--eqR, pivot_v)); //   [L,leP): < pivot_v ;    [eqR,R): == pivot_v ;

        for (int i = L; i < eqR; i++)
            if ((swap_v = data.getLong(i)) < pivot_v) data.set(i, data.set(leP++, swap_v));
            else if (swap_v == pivot_v) {
                data.set(i--, data.set(--eqR, swap_v));
            }

        for (int i = eqR,j=0; i < R; i++,j++)data.set(leP+j,data.set(i,data.getLong(leP+j)));
        eqR = leP+(R-eqR); //   [L,leP): < pivot_v ;    [leP,eqR): == pivot_v ; [eqR,R):>pivot_v
        int leftQ=qL,sameCount=R-eqR;
        while(leftQ<qR&&query_rank1[leftQ]<leP)leftQ++;
        if(qL<leftQ)getKth(L,leP,qL,leftQ);
        if(leftQ>qL&&query_rank2[leftQ-1]>=leP)leftQ--;
        while(leftQ<qR&&query_rank1[leftQ]<eqR){
            if(query_rank1[leftQ]>=leP)answer_v[leftQ]+=longToResult(pivot_v);
            if(query_rank2[leftQ]<eqR)answer_v[leftQ]+=longToResult(pivot_v);
            leftQ++;
        }
        if(leftQ>qL&&query_rank2[leftQ-1]>=eqR)leftQ--;
        if(leftQ<qR)getKth(eqR,R,leftQ,qR);
//        if (K < leP - L) return getKth(L, leP, K);
//        if (K >= (leP - L) + (R - eqR)) return getKth(leP, eqR, K - (leP - L) - (R - eqR));
    }

    public void testError(int queryN, int MULTI_QUANTILES) throws IOException {
        query_q = new double[MULTI_QUANTILES];
        answer_v = new double[MULTI_QUANTILES];
        query_rank1=new int[MULTI_QUANTILES];
        query_rank2=new int[MULTI_QUANTILES];
        final double query_delta_q = 1.0 / (MULTI_QUANTILES+1);
        for (int i = 0; i < MULTI_QUANTILES; i++){
            query_q[i] = query_delta_q * (i + 1);
            query_rank1[i] = (int) Math.floor(query_q[i] * (queryN - 1) + 1)-1; // rank1=rank2=2 for n=5&q=0.5
            query_rank2[i] = (int) Math.ceil(query_q[i] * (queryN - 1) + 1)-1;
        }

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        data = new LongArrayList(queryN);

        for (int T = 0; T < TEST_CASE; T++) {
            int L = LL[T], R = RR[T];
            data.clear();
            for(int i=L;i<R;i++)data.add(dataToLong(a[i]));
            getKth(0,queryN,0,MULTI_QUANTILES);
            for(int i=0;i<MULTI_QUANTILES;i++)answer_v[i]*=0.5;

            double[] query_a = new double[queryN];
            System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);
            for(int qid=0;qid<MULTI_QUANTILES;qid++) {
                double exact_calc_v = answer_v[qid];
                double exact_answer_v = (query_a[query_rank1[qid]]+query_a[query_rank2[qid]]) * 0.5;
//                if (DEBUG_PRINT)
//                    System.out.println("\tFINISH CALC\t\texact_calc_v=:\t" + exact_calc_v +"\tanswer_v=:\t"+exact_answer_v+ "\tqueryQ:\t"+query_q[qid]);
                    if(Math.abs(exact_answer_v-exact_calc_v)/Math.max(Math.abs(exact_answer_v),Math.abs(exact_calc_v))>1e-9)
                        System.out.println("\t\t!! CALC ERROR\n");
//                if(Math.abs(exact_quantile_v-N*query_q[qid])>1.0){
//                    System.out.println("!!!!!!!!!!ERROR!!!!");
//                    return;
//                }
            }
        }

//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("MultiQuantiles! ExactQuantile\tTEST_CASE=" + TEST_CASE + "\tDATASET:" + dataType + "\tqueryN:\t" + queryN + "\tMultiQuantiles:\t"+MULTI_QUANTILES);
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForExactQuickSelectMultiQuantiles main;

        main = new MainForExactQuickSelectMultiQuantiles();
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            main.RESULT_LINE=0;
            main.prepareA(dataType);
//            for (int queryN : new int[]{/*1000000,1500000,2000000,2500000,3000000,*/3500000,4000000,4500000,5000000/*,5500000,6000000/**/})
            for (int MULTI_QUANTILES : new int[]{/*1,1,*/1/*,4,6,10,20,40,60,100,200,400,600,1000,2000,4000,6000,10000/**/})
            for (int queryN : new int[]{100000000/*,20000000,30000000,40000000,50000000,60000000,70000000,80000000*/})
//                    for (int query_mem : new int[]{1024*128})
                    main.testError(queryN, MULTI_QUANTILES);
                main.RESULT_LINE++;
//            System.out.println("byPage & byChunkDivide\nTEST_CASE=" + TEST_CASE);
        }
//        System.out.println("\nError rate:");
//        for (String s : main.err_result)
//            System.out.println(s);
        System.out.println("\t\t\tALL_TIME:"+(new Date().getTime()-START_T));
    }
}
