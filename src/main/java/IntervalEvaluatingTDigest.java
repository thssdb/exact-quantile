import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class IntervalEvaluatingTDigest {
    int dataType;
    static int pageN = 8192;
    static int N = 55000000/pageN*pageN, pageNum=N/pageN; // CHECK IT
    public static int TEST_CASE = 1; // CHECK IT
    static double[] a;
    static TDigestForStatMerge[] workerArr;
//    String time_result = "";
    static ArrayList<String> err_result = new ArrayList<>();
    static ArrayList<String> time_result = new ArrayList<>();
    int RESULT_LINE = 0;

    Random random = new Random(233);

    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] = Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
        }
        if (dataType == 1) {
            BufferedReader reader = new BufferedReader(new FileReader(new File("1_bitcoin.csv")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN = 0;
            while ((line = reader.readLine()) != null) {
                a[cntN++] = Double.parseDouble(line);
                if (cntN == N) return;
            }
        }
        if (dataType == 2) {
            BufferedReader reader = new BufferedReader(new FileReader(new File("2_SpacecraftThruster.txt")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN = 0;
            while ((line = reader.readLine()) != null) {
                a[cntN++] = Double.parseDouble(line);
                if (cntN == N) return;
            }
        }
        if (dataType == 3) {
            BufferedReader reader = new BufferedReader(new FileReader(new File("3_taxipredition8M.txt")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN = 0;
            while ((line = reader.readLine()) != null) {
                a[cntN++] = Double.parseDouble(line);
                if (cntN == N) return;
            }
        }
        if (dataType == 4) {
            BufferedReader reader = new BufferedReader(new FileReader(new File("4_wh.csv")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN = 0;
            while ((line = reader.readLine()) != null) {
                a[cntN++] = Double.parseDouble(line);
                if (cntN == N) return;
            }
        }
    }

    public void prepareWorker(int maxSeriesByte) {
        workerArr = new TDigestForStatMerge[pageNum];
        int enoughMemByte = pageN * 10;
        for (int i = 0; i < pageNum; i++) {
            TDigestForStatMerge worker = new TDigestForStatMerge(enoughMemByte, maxSeriesByte);
            for (int j = 0; j < pageN; j++) worker.update(a[i * pageN + j]);
            worker.compactBeforeSerialization();
            workerArr[i] = worker;
            if (i == 0) {
                worker.show();
                System.out.println();
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

    public void testMergeError(int queryN, int maxMemoryByte, int maxSeriesByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time = 0, merge_time = 0;
        double err_full = 0, err_merge = 0;
        double MMP_full = 0;
        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }

        for (int T = 0; T < TEST_CASE; T++) {
            int L = LL[T], R = RR[T];
            int pageL = (L + pageN - 1) / pageN, pageR = R / pageN;
            int posL = pageL * pageN, posR = pageR * pageN;
//            System.out.println("\t\t\t"+posL+"\t"+posR);

//            merge_time-=new Date().getTime();
//            TDigestForStatMerge merge_worker = new TDigestForStatMerge(maxMemoryByte);
//            for(int i=L;i<Math.min(R,posL);i++)
//                merge_worker.update(a[i]);
//            for(int i=pageL;i<pageR;i++){
//                    merge_worker.merge(workerArr[i]);
//            }
//            for(int i=Math.max(L,posR);i<R;i++)
//                merge_worker.update(a[i]);
//            merge_time+=new Date().getTime();


            full_time -= new Date().getTime();
            TDigestForStatMerge full_worker = new TDigestForStatMerge(maxMemoryByte);
            for (int i = L; i < R; i++)
                full_worker.update(a[i]);
            full_time += new Date().getTime();

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);

            double q_add=0.0001,q_start=q_add,q_end=1-q_add,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
            for(double q=q_start;q<q_end+1e-10;q+=q_add){
                int query_rank = (int) (q * queryN);

//                double merge_v = merge_worker.quantile(q);
//                int merge_delta_rank = getValueActualRank(query_a,queryN,merge_v)-query_rank;
//                double merge_relative_err = 1.0*merge_delta_rank/(queryN);
//                err_merge+=Math.abs(merge_relative_err)/(q_count*TEST_CASE);

                double full_v = full_worker.quantile(q);

                int full_delta_rank = getDeltaRank(query_a, queryN, full_v, query_rank);
                int full_old_delta = getValueActualRank(query_a,queryN,full_v)-query_rank;//          full_delta_rank=full_old_delta;
                double full_relative_err = 1.0 * full_delta_rank / (queryN);
                err_full += Math.abs(full_relative_err) / (q_count * TEST_CASE);
                MMP_full+=1.0 *  Math.abs(Math.abs(full_old_delta)-Math.abs(full_delta_rank))  /(queryN)/(q_count*TEST_CASE);

//                int full_delta_rank = getValueActualRank(query_a, queryN, full_v) - query_rank;
//                double full_relative_err = 1.0 * full_delta_rank / (queryN);
//                err_full += Math.abs(full_relative_err) / (q_count * TEST_CASE);

//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("\t\t" + queryN +/*"\t"+err_merge+*/"\t" + err_full);//        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_mergeBuf-err_full)+"\t"+queryN*(err_merge-err_full));

        System.out.println("\t\t\t\t\t\t\tMMP DEBUG"+/*"\t"+MMP_merge+*/"\t"+MMP_full);
        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t\t" + err_full));
        time_result.set(RESULT_LINE, time_result.get(RESULT_LINE).concat("\t\t\t" + full_time));

    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START=new Date().getTime();
        IntervalEvaluatingTDigest main;

        System.out.println("interval query" + "\n");
        for (int startType=1,endType=1,dataType = startType; dataType <= endType; dataType++){ // CHECK IT
            main = new IntervalEvaluatingTDigest();
            main.prepareA(dataType);
            for (int queryN : new int[]{/*10000000, 20000000, 30000000, 40000000, */40000000})
//                for (int query_mem : new int[]{/*1024*16,1024*32,1024*64,*/1024 * 256/*,1024*256/*,1024*512,1024*1024*/})
                for(int query_mem : new int[]{1024*144,1024*168,1024*256,1024*480,1024*640,1024*1024,1024*2048}){
//                System.out.println("\n\n\t\tTDigest\tdataType:" + dataType + "\t");
                    int chunk_seri = 1 << 12; // CHECK IT
                    if (dataType == startType) {
                        err_result.add("N:" + queryN + ", " + "M:" + query_mem + ", " + "|M_c|:" + chunk_seri + "\t");
                        time_result.add("N:" + queryN + ", " + "M:" + query_mem + ", " + "|M_c|:" + chunk_seri + "\t");
                    }
//                main.prepareWorker(chunk_seri);
                    main.testMergeError(queryN, query_mem, chunk_seri); // CHECK IT
                    main.RESULT_LINE++;
                }
        }
        System.out.println("t-digest\nTEST_CASE="+TEST_CASE);
        System.out.println("\nError rate:");
        for(String s:err_result)
            System.out.println(s);
//        System.out.println("\nQuery Time:");
//        for(String s:time_result)
//            System.out.println(s);
        System.out.println("\t\tALL_TIME:"+(new Date().getTime()-START));
    }
}
