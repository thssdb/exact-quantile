import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class MainForCompareError_MERGE_FULL {
    static int pageN = 1024, pageNum=102333+12333;
    public static int TEST_CASE=128;
    static int N = pageN*pageNum;
    static double[] a;
    static KLLSketchForQuantile[] KLLArr;
    static long full_time=0,merge_time=0;
    static String time_result="";

    long seed;
    Random random = new Random(233);
    private long nextLong(){
        return random.nextLong();
    }
    private void reset(){seed=2333L;random.setSeed(233);}
    private long nextNumber(long x, int i) {
        long num=x;
        num=x&262143;
        if((i&15)<=3)
            num = nextLong();
        return num;
    }

    public static long doubleToLong(double data){
        long result = Double.doubleToLongBits(data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }
    public void prepareA(){
        a = new double[N];
        for(int i=0;i<N;i++)a[i]=random.nextGaussian();
//        for(int i=1;i<N;i++){int p = random.nextInt(i);double tmp = a[i];a[i]=a[p];a[p]=tmp;}
    }
    public void prepareKLL(int maxSeriesByte){
        KLLArr = new KLLSketchForQuantile[pageNum];
        int enoughMemByte = pageN*10;
        for(int i=0;i<pageNum;i++) {
            LongKLLSketch worker = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);
            for (int j = 0; j<pageN; j++) worker.update(doubleToLong(a[i*pageN+j]));
            worker.compactBeforeSerialization();
            KLLArr[i] = worker;
            if(i==0)worker.show();
        }
    }
    public void testErr(int queryPageNum, int maxMemoryByte, int maxSeriesByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        full_time=0;merge_time=0;
        double err_full=0,err_merge=0;
        int queryN = queryPageNum*pageN;
        double[] query_a = new double[queryN];
        for(int T=0;T<TEST_CASE;T++){
            int pageL = random.nextInt(pageNum-queryPageNum+1), pageR = pageL+queryPageNum;
            int posL = pageL*pageN, posR = pageR*pageN;
            double result_full,result_merge;

            full_time-=new Date().getTime();
            HeapLongKLLSketch full_worker = new HeapLongKLLSketch(maxMemoryByte);
            for(int i=posL;i<posR;i++)
                full_worker.update(doubleToLong(a[i]));
            full_time+=new Date().getTime();

            for(int i=posL;i<posR;i++)query_a[i-posL] = a[i];
            Arrays.sort(query_a);

            merge_time-=new Date().getTime();
            int buf_kll_num = maxMemoryByte/2/maxSeriesByte;
            List<KLLSketchForQuantile> buf_kll_list = new ArrayList<>(buf_kll_num);
            HeapLongKLLSketch merge_worker = new HeapLongKLLSketch(maxMemoryByte);
            for(int i=pageL;i<pageR;i++){
                if(buf_kll_list.size()<buf_kll_num)
                    buf_kll_list.add(KLLArr[i]);
                else{
//                    System.out.println("\t\t\t\t??!!!!!!!!????\t\t");
//                    merge_worker.show();
                    merge_worker.mergeWithTempSpace(buf_kll_list);
//                    merge_worker.show();
                    buf_kll_list.clear();
                    buf_kll_list.add(KLLArr[i]);
//                    System.out.println("\t\t\t\t??!!!!!!!!????\t\t"+Math.abs(merge_worker.getN()/2.0-merge_worker.getApproxRank(doubleToLong(merge_worker.getN()/2.0))));
                }
            }
            if(!buf_kll_list.isEmpty())
                merge_worker.mergeWithTempSpace(buf_kll_list);
            merge_time+=new Date().getTime();

//            System.out.println("\t\t\t\t??!!!!!!!!????\t\t"+Math.abs(merge_worker.getN()/2.0-merge_worker.getApproxRank(doubleToLong(merge_worker.getN()/2.0))));
//            if(T==0) {
//                full_worker.show();
//                merge_worker.show();
//                full_worker.showLevelMaxSize();
//                merge_worker.showLevelMaxSize();
//            }
            double q_start=0.05,q_end=0.96,q_add=0.02,q_count = (q_end-q_start)/q_add+1;
            for(double q=q_start;q<=q_end;q+=q_add){
                int exact_rank = (int)(q*queryN);
                double v = query_a[exact_rank];
                long full_rank = full_worker.getApproxRank(doubleToLong(v));
                long merge_rank = merge_worker.getApproxRank(doubleToLong(v));
                err_full+=Math.abs(full_rank-exact_rank)/q_count/TEST_CASE;
                err_merge+=Math.abs(merge_rank-exact_rank)/q_count/TEST_CASE;
//                System.out.println("?\t\t"+Math.abs(full_rank-exact_rank)+"\t"+Math.abs(merge_rank-exact_rank));
            }
        }
        System.out.println("\t\t"+queryPageNum+"\t"+err_merge+"\t"+err_full);
        String times=fnum.format(1.0*merge_time/TEST_CASE)+"\t"+fnum.format(1.0*full_time/TEST_CASE);
        time_result+=("\t\t\t"+times)+"\n";
    }
    public void show_time_result(){
        System.out.println(time_result);
    }
    public static void setTestCase(int tc){TEST_CASE=tc;}




    public static void main(String[] args) throws IOException{
        MainForCompareError_MERGE_FULL main;
        main = new MainForCompareError_MERGE_FULL();
        main.prepareA();
        int tmp_seri = 1<<9,tmp_mem=1<<15;
        main.prepareKLL(tmp_seri);
        for(int num=1;num<=8;num++)
            main.testErr(num,tmp_mem,tmp_seri);
        for(int num=10;num<=40960;num*=2)
            main.testErr(num,tmp_mem,tmp_seri);
        main.show_time_result();
        MainForCompareError_MERGE_FULL.setTestCase(64);
        main.testErr(81920,tmp_mem,tmp_seri);
        main.show_time_result();
        main.testErr(102333,tmp_mem,tmp_seri);
        main.show_time_result();
//        main.testErr(5,tmp_mem,tmp_seri);
//        main.testErr(10,tmp_mem,tmp_seri);
//        main.testErr(20,tmp_mem,tmp_seri);
//        main.testErr(40,tmp_mem,tmp_seri);
//        main.testErr(80,tmp_mem,tmp_seri);
//        main.testErr(200,tmp_mem,tmp_seri);
//        main.testErr(400,tmp_mem,tmp_seri);
//        main.testErr(800,tmp_mem,tmp_seri);
//        main.testErr(1600,tmp_mem,tmp_seri);
//        main.testErr(3200,tmp_mem,tmp_seri);
//        main.testErr(6400,tmp_mem,tmp_seri);
//        main.testErr(12800,tmp_mem,tmp_seri);
    }
}

// 	[8478]		[1]		[16986]		[55203]		[0]		[33811]
//	(14.29260147472178,81.0)
//	[7872]		[0]		[17326]		[13445]		[88454]
//	(12.564914867994048,76.0)