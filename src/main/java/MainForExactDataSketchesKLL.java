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
public class MainForExactDataSketchesKLL {
    int dataType=-233;
    static int startType=0,endType=4;
    static int pageN = 8192;
    static int N = pageN * 6712, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 32; // CHECK IT
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

    public int getTmpKLLK(long n,int maxMemoryByte){
        int maxItem = maxMemoryByte/8,K=0;
        for(int addK=Integer.highestOneBit(maxItem);addK>=1;addK>>=1){
            int cntK=K+addK,level=0;
            for(level=1;level<=100;level++){
                double allW=0;
                for(int i=0;i<level;i++)
                    allW+=Math.pow(2.0,i)*Math.pow(2.0/3,(level-i-1))*cntK;
                if(allW>=n)break;
            }
            double cntItem=0;
            for(int i=0;i<level;i++)
                cntItem+=Math.pow(2.0/3,(level-i-1))*cntK;
            if(cntItem<=maxItem)
                K+=cntK;
        }
//        System.out.println("\t\t\tmid iteration . KLL_K="+K);
        return K;
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

        int DataSketchesKLL_K = DataSketchesKLLForExact.calcKUnderLimitedMemory(maxMemoryByte,queryN);
//        System.out.println("\t!! DataSketchesKLL_K="+DataSketchesKLL_K);
        for (int T = 0; T < TEST_CASE; T++) {
//            System.out.println("\tTEST_ID:"+T);
//            int L = 0, R = queryN;
            int L=LL[T],R=RR[T];

            merge_page_time -= new Date().getTime();
            DataSketchesKLLForExact DataSketchesWorker = new DataSketchesKLLForExact(DataSketchesKLL_K);


            for (int i = L; i < R; i++)
                DataSketchesWorker.update(a[i]);
//                merge_page_worker.show();

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);

            double q_tmp = 1e-2, q_start = q_tmp, q_end = 1 - q_tmp;
//                q_start=0.07;q_end=0.07;
            double q_add = 1e-2, q_count = Math.floor((q_end - q_start + 1e-10) / q_add) + 1;
            double ratioPerQuery = 1.0 / (q_count * TEST_CASE);
            for (double q = q_start; q < q_end + 1e-10; q += q_add) {

                int query_rank1 = (int) Math.floor(q * queryN),query_rank2 = (int) Math.ceil(q * queryN);
                long last_n=queryN;

                double[] iterate_result = DataSketchesWorker.findResultRange(query_rank1,query_rank2);
                avg_iteration+=ratioPerQuery;
                int MMP=0;
                while(iterate_result[0]<iterate_result[1]&&iterate_result.length==2) {
                    if(++MMP>10)break;
                    avg_iteration+=ratioPerQuery;
                    DataSketchesKLLForExact cntWorker = new DataSketchesKLLForExact(getTmpKLLK(last_n,maxMemoryByte));

                    double valL=iterate_result[0],valR=iterate_result[1];
                    int CountOfLessThanValL=0;
                    for(int i=L;i<R;i++) {
                        if (a[i] >= valL && a[i] <= valR) cntWorker.update(a[i]);
                        else if (a[i] < valL) CountOfLessThanValL++;
                    }
//                    System.out.println("\t\t\t\t\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tcntSketch_N:"+cntWorker.sketch.getN());
                    int cntRank1=query_rank1-CountOfLessThanValL;
                    int cntRank2=query_rank2-CountOfLessThanValL;
//                    System.out.println("\t\t\t\t\t\tcntRank:"+cntRank1+" "+cntRank2);
                    iterate_result = cntWorker.findResultRange(cntRank1,cntRank2);
                    last_n=cntWorker.sketch.getN();
//                    System.out.println("\t\t\t\ttcntL,R:"+iterate_result[0]+","+iterate_result[1]+"\t\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tcntRank:"+cntRank1+" "+cntRank2);
                }
                double exact_quantile_v = (iterate_result[0]+iterate_result[1])*0.5;
//                System.out.println("\t\tq:"+q+"\texact_quantile="+exact_quantile_v);
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("TEST_CASE="+TEST_CASE+"\tDATASET:"+dataType+"\tqueryN:\t" + queryN+"\tmemory:\t" + maxMemoryByte+"\t\tavg_iteration:\t"+avg_iteration);
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForExactDataSketchesKLL main;

        main = new MainForExactDataSketchesKLL();
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            main.RESULT_LINE=0;
            main.prepareA(dataType);
            for (int queryN : new int[]{/*10000,100000,1000000,*/10000000})
                for (int query_mem : new int[]{1024*1,1024*2,1024*4,1024*8,1024*16,1024*32})
                        main.testError(queryN, query_mem);
                main.RESULT_LINE++;
//            System.out.println("byPage & byChunkDivide\nTEST_CASE=" + TEST_CASE);
        }
//        System.out.println("\nError rate:");
//        for (String s : main.err_result)
//            System.out.println(s);
        System.out.println("\t\t\tALL_TIME:"+(new Date().getTime()-START_T));
    }
}
