import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class IntervalEvaluatingLatencyKLL {
    int dataType = 1;
    static int pageN = 8192;
    static double[] muList = new double[]{2.0};
    static double[] sigList = new double[]{0.5,2.6,3.2};
    static int N = 55000000/pageN*pageN, pageNum=N/pageN; // CHECK IT
    public static int TEST_CASE=1; // CHECK IT
    boolean TEST_FULL = false; // CHECK IT
    static double[] a;
    static KLLSketchForQuantile[] KLLArr;
    static long[] workerMinT, workerMaxT;
    //    String time_result="";
    static ArrayList<String> err_result = new ArrayList<>();
    static ArrayList<String> time_result = new ArrayList<>();
    boolean show_time = false, show_err = true;
    static int RESULT_LINE = 0;
    Random random = new Random(233);
    static long[] aa;
    static long[] bb;
    static IntArrayList cc;

    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] = Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
        }
        BufferedReader reader = null;
        if (dataType == 1)
            reader = new BufferedReader(new FileReader(new File("1_bitcoin.csv")));
        if (dataType == 2)
            reader = new BufferedReader(new FileReader(new File("2_physiological_stress.txt")));
        if (dataType == 3)
            reader = new BufferedReader(new FileReader(new File("4_taxipredition8M.txt")));
        if (dataType == 4)
            reader = new BufferedReader(new FileReader(new File("5_wh.csv")));
        assert reader != null;
        reader.readLine(); // ignore first line.
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if (cntN == N) break;
        }
    }

    public void prepareDisorder(double mu,double sig) {
        random = new Random(233);
        aa = new long[N];
        bb = new long[N];
        cc = new IntArrayList(N);
        for (int i = 0; i < N; i++) {
            cc.add(i);
            aa[i] = i;
            bb[i] = (long) Math.round(i + Math.exp(mu + sig * random.nextGaussian()));
        }
        cc.sort((x, y) -> (Long.compare(bb[x], bb[y])));
    }

    public void prepareWorker(int maxSeriesByte) {
        KLLArr = new KLLSketchForQuantile[pageNum];
        workerMinT = new long[pageNum];
        workerMaxT = new long[pageNum];
        int enoughMemByte = pageN * 9;
        for (int i = 0; i < pageNum; i++) {
            LongKLLSketch worker = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);

            workerMinT[i] = Long.MAX_VALUE;
            workerMaxT[i] = Long.MIN_VALUE;
            for (int j = 0; j < pageN; j++) {
                int index = cc.getInt(i * pageN + j);
                worker.update(dataToLong(a[index]));
                workerMinT[i] = Math.min(workerMinT[i], index);
                workerMaxT[i] = Math.max(workerMaxT[i], index);
            }
            worker.compactBeforeSerialization();
            KLLArr[i] = worker;
//            if(i==0)worker.show();
        }
//        for(int i=0;i<100;i++)System.out.println(workerMinT[i]+"---"+workerMaxT[i]);
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


    public void testMergeError(int queryN, int maxMemoryByte, int maxSeriesByte,double sig) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time = 0, merge_time = 0;
        double err_full = 0, err_merge = 0;
        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        double  avgAvailable=0;
        for (int T = 0; T < TEST_CASE; T++) {
            int L = LL[T], R = RR[T];

            full_time -= new Date().getTime();
            HeapLongStrictKLLSketch full_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            for (int i = 0; i < N; i++)
                if(inInterval(cc.getInt(i),L,R-1))
                    full_worker.update(dataToLong(a[cc.getInt(i)]));
            full_time += new Date().getTime();

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);

            double q_add=0.0001,q_start=q_add,q_end=1-q_add,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
            for(double q=q_start;q<q_end+1e-10;q+=q_add){
                int query_rank = (int) (q * queryN);

                double full_v = longToResult(full_worker.findMinValueWithRank(query_rank));
                int full_delta_rank = getDeltaRank(query_a, queryN, full_v, query_rank);
                double full_relative_err = 1.0 * full_delta_rank / (queryN);
                err_full += Math.abs(full_relative_err) / (q_count * TEST_CASE);

//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("sig="+sig+"\t\t" + queryN + "\t" + "\t\t\t\t\t\t" + "\t" + err_full);
        //        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_mergeBuf-err_full)+"\t"+queryN*(err_merge-err_full));
        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t\t\t\t\t" + "\t" + err_full));
//        time_result.set(RESULT_LINE, time_result.get(RESULT_LINE).concat("\t\t\t\t\t\t" + merge_time+(TEST_FULL?("\t"+full_time):"")));

    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        IntervalEvaluatingLatencyKLL main;

        System.out.println("UnSeqData. OnlineKLL. interval query" + "\n");
        for (int dataType = 1; dataType <= 1; dataType++) { // CHECK IT
            main = new IntervalEvaluatingLatencyKLL();
            main.prepareA(dataType);
            for (double mu : muList)
                for (double sig : sigList) {
                main.prepareDisorder(mu,sig);
                err_result.add("mu:" + mu + ", " + "sig:" + sig + "\t");
                for (int queryN : new int[]{40000000})
                    for (int query_mem : new int[]{1024 * 256})
                        for (int chunk_seri : new int[]{/*128,256,512,*/4096/*,128,256,512,1024/*,2048,4096,8192*/})
                            main.testMergeError(queryN, query_mem, chunk_seri,sig);
                IntervalEvaluatingLatencyKLL.RESULT_LINE++;
            }
            System.out.println("Online KLL\nTEST_CASE=" + TEST_CASE);
            System.out.println("\nError rate:");
            for (String s : err_result)
                System.out.println(s);
        }
    }
}
