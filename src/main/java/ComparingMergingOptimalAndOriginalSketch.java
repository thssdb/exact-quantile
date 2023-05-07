import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class ComparingMergingOptimalAndOriginalSketch {
    int dataType = 1;
    static int pageN = 8192*16;
    static int N = 81000000/pageN*pageN+pageN, pageNum=N/pageN; // CHECK IT
    public static int TEST_CASE=128; // CHECK IT
    static double[] a;
    static KLLSketchForQuantile[] OptimalArr;
    static HeapLongStrictKLLSketch[] OriginalArr;
    String time_result="";
    static ArrayList<String> err_result = new ArrayList<>();
    boolean show_time = false,show_err=true;
    int RESULT_LINE=0;

    Random random = new Random(233);

    public void prepareA(int dataType) throws IOException {
        if(a==null)a = new double[N];
        this.dataType = dataType;

        if(dataType==0){
            for(int i=0;i<N;i++)a[i]=
                random.nextGaussian();
                //Math.pow(-1,random.nextInt(2))*Math.pow(10.0,(2*Math.pow(random.nextDouble(),2)-1)*300);
        }
        if(dataType==1){
            BufferedReader reader = new BufferedReader(new FileReader(new File("1_bitcoin.csv")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN= 0;
            while((line=reader.readLine())!=null){
                a[cntN++] = Double.parseDouble(line);
                if(cntN==N)return;
            }
        }
        if(dataType==2){
            BufferedReader reader = new BufferedReader(new FileReader(new File("2_SpacecraftThruster.txt")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN= 0;
            while((line=reader.readLine())!=null){
                a[cntN++] = Double.parseDouble(line);
                if(cntN==N)return;
            }
        }
        if(dataType==3){
            BufferedReader reader = new BufferedReader(new FileReader(new File("3_taxipredition8M.txt")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN= 0;
            while((line=reader.readLine())!=null){
                a[cntN++] = Double.parseDouble(line);
                if(cntN==N)return;
            }
        }
        if(dataType==4){
            BufferedReader reader = new BufferedReader(new FileReader(new File("4_wh.csv")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN= 0;
            while((line=reader.readLine())!=null){
                a[cntN++] = Double.parseDouble(line);
                if(cntN==N)return;
            }
        }
    }
    public void prepareWorker(int maxSeriesByte){
        long ST;
        OptimalArr = new KLLSketchForQuantile[pageNum];
        OriginalArr = new HeapLongStrictKLLSketch[pageNum];

        ST = new Date().getTime();
        int enoughMemByte = pageN*8;
        for(int i=0;i<pageNum;i++) {
            LongKLLSketch worker = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);
            for (int j = 0; j<pageN; j++) worker.update(dataToLong(a[i*pageN+j]));
            worker.compactBeforeSerialization();
            OptimalArr[i] = worker;
            if(i==0)worker.show();
        }
        if(show_time) {
        System.out.print("\t\t"+(new Date().getTime()-ST));
        }

        ST = new Date().getTime();
        for(int i=0;i<pageNum;i++) {
            HeapLongStrictKLLSketch original_worker = new HeapLongStrictKLLSketch(maxSeriesByte);
            for (int j = 0; j<pageN; j++) original_worker.update(dataToLong(a[i*pageN+j]));
            OriginalArr[i] = original_worker;
            if(i==0)original_worker.show();
        }
        if(show_time) {
            System.out.print("\t" + (new Date().getTime() - ST) + "\n");
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


    private long dataToLong(double data)  {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }

    private double longToResult(long result)  {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }

    public void testMergeError(int queryN, int maxMemoryByte, int maxSeriesByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long original_time=0, optimal_time=0;
        double err_original=0,err_optimal=0;
        double[] query_a = new double[queryN];


        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = 0;//random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }

        for(int T=0;T<TEST_CASE;T++){
            prepareWorker(maxSeriesByte);
            int L = LL[T], R = RR[T];
            int pageL = (L+pageN-1)/pageN, pageR = R/pageN;
            int posL = pageL*pageN, posR = pageR*pageN;
//            System.out.println("\t\t\t"+posL+"\t"+posR);

            optimal_time-=new Date().getTime();
            int buf_kll_num = 1;
            HeapLongStrictKLLSketch optimal_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            for(int i=L;i<Math.min(R,posL);i++)
                optimal_worker.update(dataToLong(a[i]));
            for(int i=pageL;i<pageR;i++)
                optimal_worker.mergeWithTempSpace(OptimalArr[i]);
            for(int i=Math.max(L,posR);i<R;i++)
                optimal_worker.update(dataToLong(a[i]));
            optimal_time+=new Date().getTime();



            original_time-=new Date().getTime();
            HeapLongStrictKLLSketch original_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            for(int i=L;i<Math.min(R,posL);i++)
                original_worker.update(dataToLong(a[i]));
            for(int i=pageL;i<pageR;i++)
                original_worker.mergeWithTempSpace(OriginalArr[i]);
            for(int i=Math.max(L,posR);i<R;i++)
                original_worker.update(dataToLong(a[i]));
            original_time+=new Date().getTime();

            optimal_worker.show();
            original_worker.show();

            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R-L);
            Arrays.sort(query_a);

            double q_start=0.0001,q_end=1-q_start,q_add=0.0001,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
            for(double q=q_start;q<q_end+1e-10;q+=q_add){
                int query_rank = (int)(q*queryN);

                double optimal_v = longToResult(optimal_worker.findMinValueWithRank(query_rank));
                int optimal_delta_rank = getValueActualRank(query_a,queryN,optimal_v)-query_rank;
                double optimal_relative_err = 1.0*optimal_delta_rank/(queryN);
                err_optimal+=Math.abs(optimal_relative_err)/(q_count*TEST_CASE);

                double original_v = longToResult(original_worker.findMinValueWithRank(query_rank));
                int original_delta_rank = getValueActualRank(query_a,queryN,original_v)-query_rank;
                double original_relative_err = 1.0*original_delta_rank/(queryN);
                err_original+=Math.abs(original_relative_err)/(q_count*TEST_CASE);

//                System.out.println("?\t\toriginal:"+original_v+" delta:"+original_delta_rank+"\t\toptimal:"+optimal_v+" delta:"+optimal_delta_rank);
            }
        }
        if(show_err) {
            System.out.print("\t\t\t" + err_optimal + "\t" + err_original);
            err_result.set(RESULT_LINE,err_result.get(RESULT_LINE).concat("\t\t\t"+err_optimal+"\t"+err_original));
        }
        if(TEST_CASE!=0&&show_time) {
            System.out.print("\t" + optimal_time / TEST_CASE);
            System.out.print("\t" + original_time / TEST_CASE);
            System.out.println();
        }
        //        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_optimalBuf-err_original)+"\t"+queryN*(err_optimal-err_original));

    }
    public void show_time_result(){
        System.out.println(time_result);
    }
    public static void setTestCase(int tc){TEST_CASE=tc;}




    public static void main(String[] args) throws IOException{
        long START=new Date().getTime();
        ComparingMergingOptimalAndOriginalSketch main;
//        main = new MainForMergeStatErrorKLL();
//        main.prepareA();
//        int tmp_seri = 1<<9,tmp_mem=1<<15;
//        main.prepareKLL(tmp_seri);
////        for(int num=1;num<=8;num++)
////            main.testKLL(num,tmp_mem,tmp_seri);
//        for(int num=1;num<=pageNum;num*=2)
//            main.testKLL(num,tmp_mem,tmp_seri);
//        main.show_time_result();

        System.out.println("interval query"+"\n");
        for (int startType=2,endType=2,dataType = startType; dataType <= endType; dataType++){ // CHECK IT
            main = new ComparingMergingOptimalAndOriginalSketch();
            main.prepareA(dataType);
            for(int queryN : new int[]{/*10000000,20000000,30000000,40000000,*/50000000})
                for(int query_mem : new int[]{/*1024*16,1024*32,1024*64,1024*128/*,*/1024*256*8/*,1024*512,1024*1024*/})
                    for(int chunk_seri : new int[]{/*128,256,512,1024/*,128,256,512,1024/*,2048,4096,8192*//*256*16,512*16,1024*16,2048*16,4096*16*/1024*2,1024*4,1024*8,1024*16,1024*24,1024*32,1024*40,1024*48,1024*56,1024*64/*,1024*24,1024*32,1024*40,1024*48,1024*56,1024*64*/}) {
//                        System.out.println("\n\n\t\tKLL\tdataType:" + dataType + "\t");
                        int tmp_seri = chunk_seri, tmp_mem = query_mem; // CHECK IT
//                        System.out.println("\tpageN:" + pageN + "\t|summary|:" + tmp_seri + "\t|Memory|:" + tmp_mem);
//                        main.prepareWorker(tmp_seri);
                        if(dataType==startType)
                            err_result.add("N:"+queryN+", "+"M:\t"+query_mem+", "+"|M_c|:\t"+chunk_seri+"\t");
                        main.testMergeError(queryN, tmp_mem, tmp_seri); // CHECK IT
                        System.out.println("");
                        main.RESULT_LINE++;
//                        main.show_time_result();
                }
            }
        for(String s:err_result)
            System.out.println(s);
        System.out.println("\t\tALL_TIME:"+(new Date().getTime()-START));
    }
}
