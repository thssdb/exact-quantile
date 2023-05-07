import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;

public class IntervalEvaluatingMergingSampling {
    int dataType;
    static int N = 50000000; // CHECK IT
    static int pageN = 1 << 12, pageNum = N / pageN;
    public static int TEST_CASE = 4; // CHECK IT
    static double[] a;
    static SamplingHeapForStatMerge[] workerArr;
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
            BufferedReader reader = new BufferedReader(new FileReader(new File("2_physiological_stress.txt")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN = 0;
            while ((line = reader.readLine()) != null) {
                a[cntN++] = Double.parseDouble(line);
                if (cntN == N) return;
            }
        }
        if (dataType == 3) {
            BufferedReader reader = new BufferedReader(new FileReader(new File("4_taxipredition8M.txt")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN = 0;
            while ((line = reader.readLine()) != null) {
                a[cntN++] = Double.parseDouble(line);
                if (cntN == N) return;
            }
        }
        if (dataType == 4) {
            BufferedReader reader = new BufferedReader(new FileReader(new File("5_wh.csv")));
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
        workerArr = new SamplingHeapForStatMerge[pageNum];
        int enoughMemByte = pageN * 10;
        for (int i = 0; i < pageNum; i++) {
            SamplingHeapForStatMerge worker = new SamplingHeapForStatMerge(enoughMemByte, maxSeriesByte);
            for (int j = 0; j < pageN; j++) worker.update(dataToLong(a[i * pageN + j]));
            worker.compactBeforeSerialization();
            workerArr[i] = worker;
            if (i == 0) {
                worker.show();
                System.out.println();
            }
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


    private long dataToLong(double data) {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }

    private double longToResult(long result) {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }

    public void testMergeError(int L, int R, int maxMemoryByte, int maxSeriesByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time = 0, merge_time = 0;
        double err_full = 0, err_merge = 0;
        int queryN = R - L;
        double[] query_a = new double[queryN];

        for (int T = 0; T < TEST_CASE; T++) {
            int pageL = (L + pageN - 1) / pageN, pageR = R / pageN;
            int posL = pageL * pageN, posR = pageR * pageN;
//            System.out.println("\t\t\t"+posL+"\t"+posR);

//            merge_time-=new Date().getTime();
//            SamplingHeapForStatMerge merge_worker = new SamplingHeapForStatMerge(maxMemoryByte);
//            for(int i=L;i<Math.min(R,posL);i++)
//                merge_worker.update(dataToLong(a[i]));
//            for(int i=pageL;i<pageR;i++){
//                    merge_worker.merge(workerArr[i]);
//            }
//            for(int i=Math.max(L,posR);i<R;i++)
//                merge_worker.update(dataToLong(a[i]));
//            merge_time+=new Date().getTime();


            full_time -= new Date().getTime();
            SamplingHeapForStatMerge full_worker = new SamplingHeapForStatMerge(maxMemoryByte);
            for (int i = L; i < R; i++)
                full_worker.update(dataToLong(a[i]));
            full_time += new Date().getTime();

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);

            double q_start = 0.01, q_end = 0.99, q_add = 0.005, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
                int query_rank = (int) (q * queryN);

//                double merge_v = longToResult(merge_worker.quantile(q));
//                int merge_delta_rank = getValueActualRank(query_a,queryN,merge_v)-query_rank;
//                double merge_relative_err = 1.0*merge_delta_rank/(queryN);
//                err_merge+=Math.abs(merge_relative_err)/(q_count*TEST_CASE);

                double full_v = longToResult(full_worker.quantile(q));
                int full_delta_rank = getValueActualRank(query_a, queryN, full_v) - query_rank;
                double full_relative_err = 1.0 * full_delta_rank / (queryN);
                err_full += Math.abs(full_relative_err) / (q_count * TEST_CASE);

//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("\t\t" + queryN +/*"\t"+err_merge+*/"\t" + err_full);//        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_mergeBuf-err_full)+"\t"+queryN*(err_merge-err_full));
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
        IntervalEvaluatingMergingSampling main;

        System.out.println("interval query" + "\n");
        for (int dataType = 1; dataType <= 4; dataType++) { // CHECK IT
            main = new IntervalEvaluatingMergingSampling();
            main.prepareA(dataType);
            for (int i : new int[]{10000000,20000000,30000000,40000000,50000000/**/})
                for (int query_mem : new int[]{/*1024 * 16, 1024 * 32, 1024 * 64, */1024 * 128/*, 1024 * 256/*,1024*512,1024*1024*/}) {
//                    System.out.println("\n\n\t\tTDigest\tdataType:" + dataType + "\t");
                    int chunk_seri = 1 << 10; // CHECK IT
//                System.out.println("\tpageN:"+pageN+"\t|summary|:"+tmp_seri+"\t|Memory|:"+tmp_mem);
//                main.prepareWorker(tmp_seri);
                    if (dataType == 1) {
                        err_result.add("N:" + i + ", " + "M:" + query_mem + ", " + "|M_c|:" + chunk_seri + "\t");
                        time_result.add("N:" + i + ", " + "M:" + query_mem + ", " + "|M_c|:" + chunk_seri + "\t");
                    }
                    main.testMergeError(0, i, query_mem, chunk_seri); // CHECK IT
//                    main.show_time_result();
                    main.RESULT_LINE++;
                }
        }
        System.out.println("random sampling\nTEST_CASE="+TEST_CASE);
        System.out.println("\nError rate:");
        for(String s:err_result)
            System.out.println(s);
        System.out.println("\nQuery Time:");
        for(String s:time_result)
            System.out.println(s);
    }
}