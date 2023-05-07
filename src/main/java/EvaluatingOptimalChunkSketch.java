import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class EvaluatingOptimalChunkSketch {
    static boolean errRecord = false;
    static int dataType = 1;
    static int pageN = 8192;
    static int N = 55000000/pageN*pageN, pageNum=N/pageN; // CHECK IT
    public static int TEST_CASE = 1024; // CHECK IT
    static double[] a;
    static KLLSketchForQuantile[] KLLArr;
    static ArrayList<String> err_result = new ArrayList<>();
    int RESULT_LINE = 0;
    String time_result = "";

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
        KLLArr = new KLLSketchForQuantile[pageNum];
        int enoughMemByte = pageN * 10;
        for (int i = 0; i < pageNum; i++) {
            LongKLLSketch worker = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);
            for (int j = 0; j < pageN; j++) worker.update(dataToLong(a[i * pageN + j]));
            worker.compactBeforeSerialization();
            KLLArr[i] = worker;
            if (i == 0) worker.show();
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

    public void testSeriError(int maxSeriesByte) throws IOException {
        DecimalFormat fnum = new DecimalFormat("#0.00");
        double err_optimal = 0, err_stream = 0;
        int queryN = pageN;
        double[] query_a = new double[queryN];

        for (int T = 0; T < TEST_CASE; T++) {
            int L = random.nextInt(N - pageN + 1), R = L + pageN;

            int buf_kll_num = 1;

            LongKLLSketch optimal_worker = new LongKLLSketch(pageN, pageN * 10, maxSeriesByte);
            for (int i = L; i < R; i++)
                optimal_worker.update(dataToLong(a[i]));
            optimal_worker.compactBeforeSerialization();


            HeapLongStrictKLLSketch stream_worker = new HeapLongStrictKLLSketch(maxSeriesByte);
            for (int i = L; i < R; i++)
                stream_worker.update(dataToLong(a[i]));

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);

//            optimal_worker.showNum();
//            stream_worker.showNum();

            double q_add=0.0001,q_start=q_add,q_end=1-q_add,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
            for(double q=q_start;q<q_end+1e-10;q+=q_add){
                int query_rank = (int) (q * queryN);

                double optimal_v = longToResult(optimal_worker.findMinValueWithRank(query_rank));
                int optimal_delta_rank = getDeltaRank(query_a, queryN, optimal_v, query_rank);
                double optimal_relative_err = 1.0 * optimal_delta_rank / (queryN);
                err_optimal += Math.abs(optimal_relative_err) / (q_count * TEST_CASE);

                double stream_v = longToResult(stream_worker.findMinValueWithRank(query_rank));
                int stream_delta_rank = getDeltaRank(query_a, queryN, stream_v, query_rank);
                double stream_relative_err = 1.0 * stream_delta_rank / (queryN);
                err_stream += Math.abs(stream_relative_err) / (q_count * TEST_CASE);

//                System.out.println("?\t\toptimal:"+optimal_v+" delta:"+optimal_delta_rank+"\t\tstream:"+stream_v+" delta:"+stream_delta_rank);
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("\t\t\t" + maxSeriesByte + "\t" + err_optimal + "\t" + err_stream);//        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_mergeBuf-err_full)+"\t"+queryN*(err_merge-err_full));

        err_result.set(RESULT_LINE,err_result.get(RESULT_LINE).concat("\t\t\t"+err_optimal+"\t"+err_stream));
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START=new Date().getTime();
        EvaluatingOptimalChunkSketch main;
//        main = new MainForMergeStatErrorKLL();
//        main.prepareA();
//        int tmp_seri = 1<<9,tmp_mem=1<<15;
//        main.prepareKLL(tmp_seri);
////        for(int num=1;num<=8;num++)
////            main.testKLL(num,tmp_mem,tmp_seri);
//        for(int num=1;num<=pageNum;num*=2)
//            main.testKLL(num,tmp_mem,tmp_seri);
//        main.show_time_result();

        System.out.println("interval query" + "\n");
        for (int dataType = 1; dataType <= 4; dataType++) {
            main = new EvaluatingOptimalChunkSketch();
            main.prepareA(dataType);

            System.out.println("\n\n\t\tKLL\tdataType:" + dataType + "\t");
            for (int chunk_seri : new int[]{512,1024,2048,4096,8192}) {
                if (dataType == 1) {
                    err_result.add("PageN:" + pageN + "|M_c|:" + chunk_seri + "\t");
                }
//                    System.out.println("\tpageN:" + pageN + "\t|summary|:" + tmp_seri);
                main.testSeriError(chunk_seri); // CHECK IT
                main.RESULT_LINE++;
//                    main.show_time_result();
            }
        }
        System.out.println("Comparing Chunk Sketch: optimal & normal\nTEST_CASE="+TEST_CASE);
        System.out.println("\nError rate for single chunk data:");
        for(String s:err_result)
            System.out.println(s);
        System.out.println("\t\tALL_TIME:"+(new Date().getTime()-START));
    }
}
