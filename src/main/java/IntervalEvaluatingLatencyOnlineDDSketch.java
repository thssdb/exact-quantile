import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntDoubleImmutablePair;
import it.unimi.dsi.fastutil.ints.IntDoublePair;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class IntervalEvaluatingLatencyOnlineDDSketch {
    int dataType = 1;
    static int LSM_T=30;
    static int pageN = 8192;
    static double[] muList = new double[]{3.0};
    static double[] sigList = new double[]{0.0,1.0,2.1,2.4,2.6};
    static int N = 55000000/pageN*pageN, pageNum=N/pageN; // CHECK IT
    public static int TEST_CASE=1; // CHECK IT
    boolean TEST_FULL = false; // CHECK IT
    static double[] a;
    static double DDSketch_ALPHA;
    static KLLSketchForQuantile[] KLLArr;
    static long[] workerMinT, workerMaxT;
    //    String time_result="";
    static ArrayList<String> err_result = new ArrayList<>();
    static ArrayList<String> time_result = new ArrayList<>();
    boolean show_time = false, show_err = true;
    int RESULT_LINE = 0;
    Random random = new Random(233);
    static long[] aa;
    static long[] bb;
    static IntArrayList cc;
    ObjectArrayList<IntDoublePair>allCompactedData;

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

    public void prepareWorker() {

        int sketchNum = 1,SSTlv=0;
        while(sketchNum*LSM_T<=pageNum) {
            sketchNum*=LSM_T;
            SSTlv++;
        }
        allCompactedData = new ObjectArrayList<>();
//        System.out.println("\t\t??\t"+lsmNode.size());
        int cntCompactedPos=0;
        for(int lsmLV=SSTlv;lsmLV>=0;lsmLV--){
            int chunkNum = (int)Math.pow(LSM_T,lsmLV);
            while(cntCompactedPos+chunkNum<=pageNum){
                ObjectArrayList<IntDoublePair>compactedData =new ObjectArrayList<>();
                for(int i=0;i<chunkNum*pageN;i++){
                    int pos =cntCompactedPos*pageN+i;
                    IntDoublePair pair = new IntDoubleImmutablePair(cc.getInt(pos),a[cc.getInt(pos)]);
                    compactedData.add(pair);
//                    if(compactedData.size()%100==0)System.out.println("\t\t\t\tdata:"+pair.firstInt()+"\t\tv="+pair.secondDouble());
                }
//                System.out.println("\t\t\t\t\tre-sorted data   range:\t"+cntCompactedPos*pageN+"......"+(cntCompactedPos*pageN+chunkNum*pageN));
                compactedData.sort(Comparator.comparingInt(IntDoublePair::firstInt));
                allCompactedData.addAll(compactedData);
//                System.out.println("\t\t\t\t\t\tsst mnTmxT:"+compactedData.get(0).firstInt()+"..."+compactedData.get(compactedData.size()-1).firstInt());

                cntCompactedPos+=chunkNum;
            }
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

            DDSketchForQuantile full_worker = new DDSketchForQuantile(DDSketch_ALPHA,maxMemoryByte/42);
            for (int i = L; i < R; i++)
                if(inInterval(allCompactedData.get(i).firstInt(),L,R-1)) {
                    full_worker.insert(allCompactedData.get(i).secondDouble());
                }
            full_time += new Date().getTime();

//            System.out.println("\t\t\tquery l,r:"+L+"\t\t"+R);
//            full_worker.show();

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);

            double q_add=0.0001,q_start=q_add,q_end=1-q_add,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
            for(double q=q_start;q<q_end+1e-10;q+=q_add){
                int query_rank = (int) (q * queryN);

                double full_v = full_worker.getQuantile(q);
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
        long START=new Date().getTime();
        IntervalEvaluatingLatencyOnlineDDSketch main;

        System.out.println("UnSeqData. DDSketch. interval query" + "\n");
        for (int startType=1,endType=4,dataType = startType; dataType <= endType; dataType++){ // CHECK IT
            main = new IntervalEvaluatingLatencyOnlineDDSketch();
            main.prepareA(dataType);
            for (double mu : muList)
                for (double sig : sigList) {
                main.prepareDisorder(mu,sig);
                main.prepareWorker();
                for (int queryN : new int[]{40000000})
                    for (int query_mem : new int[]{1024 * 256})
                        for (int chunk_seri : new int[]{/*128,256,512,*/4096/*,128,256,512,1024/*,2048,4096,8192*/}) {

                            int dataset_V = 100000,limit = query_mem/42;
                            DDSketch_ALPHA = Math.pow(10,Math.log10(dataset_V)/limit)-1;
                            if(dataType==startType) {
                                err_result.add("sig:\t"+sig+"\tN:" + queryN + ", " +  ", " + "M:" + query_mem / 1024 + "\t");
                                time_result.add("N:" + queryN + ", " + "M:" + query_mem / 1024 + "\t");
                            }
                            main.testMergeError(queryN, query_mem, chunk_seri, sig);
                            main.RESULT_LINE++;
                        }
            }
        }
        System.out.println("DDSketch\nTEST_CASE=" + TEST_CASE);
        System.out.println("\nError rate:");
        for (String s : err_result)
            System.out.println(s);
        System.out.println("\t\tALL_TIME:"+(new Date().getTime()-START));
    }
}
