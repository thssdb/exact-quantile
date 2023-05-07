import it.unimi.dsi.fastutil.doubles.DoubleArrayList;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.util.*;

public class MainForMergeStatErrorTDigest {
    static boolean errRecord = true;
    static int dataType = 0;
    static int pageN = 1<<10, pageNum=32768*2+5000;
    public static int TEST_CASE=24;
    static int N = pageN*pageNum;
    static double[] a;
    static TDigestForStatMerge[] workerArr;
    static long full_time=0,merge_time=0;
    String time_result="";

    Random random = new Random(233);

    public void prepareA(){
        a = new double[N];
//        for(int i=0;i<N;i++)a[i]=i;for(int i=1;i<N;i++){int p = random.nextInt(i);double tmp = a[i];a[i]=a[p];a[p]=tmp;}
        if(dataType==0){
            for(int i=0;i<N;i++)a[i]=random.nextGaussian();
        }
        if(dataType==1){
            for(int i=0;i<N;i++)a[i]=random.nextGaussian();
            Arrays.sort(a,0,N);
        }
        if(dataType==2){
            for(int i=0;i<N;i++)a[i]=Math.pow(-1,random.nextInt(2))*Math.pow(10.0,(2*Math.pow(random.nextDouble(),2)-1)*300);
        }
    }
    public void prepareWorker(int maxSeriesByte){
        workerArr = new TDigestForStatMerge[pageNum];
        int enoughMemByte = pageN*10;
        for(int i=0;i<pageNum;i++) {
            TDigestForStatMerge worker = new TDigestForStatMerge(enoughMemByte, maxSeriesByte);
            for (int j = 0; j<pageN; j++) worker.update(a[i*pageN+j]);
            worker.compactBeforeSerialization();
            workerArr[i] = worker;
            if(i==0){worker.show();System.out.println();}
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

    public void testMergeError(int queryPageNum, int maxMemoryByte, int maxSeriesByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        full_time=0;merge_time=0;
        double err_full=0,err_merge=0;
        int queryN = queryPageNum*pageN;
        double[] query_a = new double[queryN];
        DoubleArrayList err_record_merge = new DoubleArrayList();
        DoubleArrayList err_record_full = new DoubleArrayList();

        for(int T=0;T<TEST_CASE;T++){
            int pageL = random.nextInt(pageNum-queryPageNum+1), pageR = pageL+queryPageNum;
            int posL = pageL*pageN, posR = pageR*pageN;
            double result_full,result_merge;


            merge_time-=new Date().getTime();
            int buf_kll_num = 1;//maxMemoryByte/maxSeriesByte/4;
            TDigestForStatMerge merge_worker = new TDigestForStatMerge(maxMemoryByte);
            for(int i=pageL;i<pageR;i++){
                merge_worker.merge(workerArr[i]);
            }
            merge_time+=new Date().getTime();


            full_time-=new Date().getTime();
            TDigestForStatMerge full_worker = new TDigestForStatMerge(maxMemoryByte);
            for(int i=posL;i<posR;i++)
                full_worker.update(a[i]);
            full_time+=new Date().getTime();


//            System.out.println("\t\t\t\t??!!!!!!!!????\t\t"+Math.abs(merge_worker.getN()/2.0-merge_worker.getApproxRank(doubleToLong(merge_worker.getN()/2.0))));
            if(T==0) {
//                merge_worker.show();
//                full_worker.show();
//                full_worker.showLevelMaxSize();
//                merge_worker.showLevelMaxSize();
            }
            if (posR - posL >= 0) System.arraycopy(a, posL, query_a, 0, posR - posL);
            Arrays.sort(query_a);
            double q_start=0.01,q_end=0.99,q_add=0.01,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
            for(double q=q_start;q<q_end+1e-10;q+=q_add){

                int query_rank = (int)(q*queryN);
                double merge_v = merge_worker.quantile(q);
                int merge_delta_rank = getValueActualRank(query_a,queryN,merge_v)-query_rank;
                double merge_relative_err = 1.0*merge_delta_rank/(queryPageNum*pageN);
                err_merge+=Math.abs(merge_relative_err)/(q_count*TEST_CASE);
                err_record_merge.add(merge_relative_err);

                double full_v = full_worker.quantile(q);
                int full_delta_rank = getValueActualRank(query_a,queryN,full_v)-query_rank;
                double full_relative_err = 1.0*full_delta_rank/(queryPageNum*pageN);
                err_full+=Math.abs(full_relative_err)/(q_count*TEST_CASE);
                err_record_full.add(full_relative_err);

//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
            }
        }
        System.out.println("\t\t"+queryPageNum+"\t"+err_merge+"\t"+err_full);
        String times=fnum.format(1.0*merge_time/TEST_CASE)+"\t"+fnum.format(1.0*full_time/TEST_CASE);
        time_result+=("\t\t"+queryPageNum+"\t"+times)+"\n";
        if(errRecord) {
            err_record_merge.sort(Double::compare);
            err_record_full.sort(Double::compare);
            File f1 = new File("./" + queryPageNum + "_" + maxMemoryByte + "_" + maxSeriesByte + "_"+dataType+"_TDigest");
            FileOutputStream fos1 = new FileOutputStream(f1);
            OutputStreamWriter w1 = new OutputStreamWriter(fos1, StandardCharsets.UTF_8);
            for (int i = 0; i < err_record_merge.size(); i++)
                w1.append("\t").append(String.valueOf(err_record_merge.getDouble(i))).append("\t").append(String.valueOf(err_record_full.getDouble(i))).append("\n");
            w1.close();
            fos1.close();
        }
    }
    public void show_time_result(){
        System.out.println("\tTimeResult:");
        System.out.println(time_result);
    }
    public void setTestCase(int tc){TEST_CASE=tc;}


    public void checkTDigestAvg(){
        long START_TIME = new Date().getTime();
        int maxMemoryByte = 1<<20,maxSeriByte = 1<<9;
//        TDigestForStatMerge worker = new TDigestForStatMerge(maxMemoryByte,maxSeriByte);
//        for(int i=0;i<9200;i++)worker.update(i);
//        worker.show();
//        worker.compactBeforeSerialization();
//        worker.show();
        for(int T=0;T<1;T++) {
            TDigestForStatMerge worker_avg = new TDigestForStatMerge(maxMemoryByte, maxSeriByte);
            for (int i = 0; i < pageN; i++)
                worker_avg.update(i);
//        worker_avg.sortCluster();
//        worker_avg.show();
            worker_avg.compactBeforeSerialization();
            worker_avg.show();
        }
        System.out.println("\t\ttime:"+(new Date().getTime()-START_TIME));
    }
    public void checkTDigestMid(){
        long START_TIME = new Date().getTime();
        int maxMemoryByte = 1<<10,maxSeriByte = 1<<10;
        TDigestMidForStatMerge worker_mid = new TDigestMidForStatMerge(maxMemoryByte,maxSeriByte);
        for(int i=0;i<(1<<20);i++)
            worker_mid.update(i);
        worker_mid.sortCluster();
        worker_mid.show();
        System.out.println("\t\ttime:"+(new Date().getTime()-START_TIME));
    }
    public void checkTDigestSpeed(){
        int maxMemoryByte = 1<<15,maxSeriByte = 1<<10;
        long START_TIME;


        START_TIME = new Date().getTime();
        TDigestTedForStatMerge worker_ted = new TDigestTedForStatMerge(maxMemoryByte,maxSeriByte);
        for(int i=0;i<(1<<24);i++)
            worker_ted.update(i);
        System.out.println("\t\ttime:"+(new Date().getTime()-START_TIME));
        START_TIME = new Date().getTime();
        worker_ted = new TDigestTedForStatMerge(maxMemoryByte,maxSeriByte);
        for(int i=0;i<(1<<24);i++)
            worker_ted.update(i);
        System.out.println("\t\ttime:"+(new Date().getTime()-START_TIME));
    }

    public static void main(String[] args) throws IOException{
        MainForMergeStatErrorTDigest main;
        for(int dataType=0;dataType<=2;dataType++) {
            MainForMergeStatErrorTDigest.dataType = dataType;
            System.out.println("\n\n\t\tT-Digest\tdataType:"+dataType+"\t");
            main = new MainForMergeStatErrorTDigest();
            main.prepareA();
            int tmp_seri = 1 << 9, tmp_mem = 1 << 15;
            main.prepareWorker(tmp_seri);
            for (int num = 1; num <= pageNum; num *= 2)
                main.testMergeError(num, tmp_mem, tmp_seri);
            main.show_time_result();
        }
//        main.checkTDigestAvg();
//        main.checkTDigestSpeed();
    }
}
