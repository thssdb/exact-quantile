import it.unimi.dsi.fastutil.doubles.DoubleArrayList;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.util.*;

public class MainForMergeStatErrorKLL {
    static boolean errRecord = false;
    static int dataType = 0;
    static int pageN = 1<<12, pageNum=/*32768*2+5000*/20000;
    public static int TEST_CASE=64;
    static int N = pageN*pageNum;
    static double[] a;
    static KLLSketchForQuantile[] KLLArr;
    static final double memoryForMergeBuffer=0.25;
    String time_result="";

    Random random = new Random(233);

    public void prepareA(){
        a = new double[N];
//        for(int i=0;i<N;i++)a[i]=i;for(int i=1;i<N;i++){int p = random.nextInt(i);double tmp = a[i];a[i]=a[p];a[p]=tmp;}
//        if(dataType==0){
//            for(int i=0;i<N;i++)a[i]=random.nextGaussian();
//        }
//        if(dataType==1){
//            for(int i=0;i<N;i++)a[i]=random.nextGaussian();
//            Arrays.sort(a,0,N);
//        }
        if(dataType==0){
            for(int i=0;i<N;i++)a[i]=Math.pow(-1,random.nextInt(2))*Math.pow(10.0,(2*Math.pow(random.nextDouble(),2)-1)*300);
        }
    }
    public void prepareWorker(int maxSeriesByte){
        KLLArr = new KLLSketchForQuantile[pageNum];
        int enoughMemByte = pageN*10;
        for(int i=0;i<pageNum;i++) {
            LongKLLSketch worker = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);
            for (int j = 0; j<pageN; j++) worker.update(dataToLong(a[i*pageN+j]));
            worker.compactBeforeSerialization();
            KLLArr[i] = worker;
            if(i==0)worker.show();
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

    private long approximatePointAvgError(int pageKLLNum, int maxMemoryByte, int SKETCH_BYTE) {
        long TOT_N = (long) pageKLLNum * pageN;
        int MemNum = maxMemoryByte/8;
        return (long) Math.ceil(4/3.0 * TOT_N / MemNum);
    }
    private long approximateMergeBufAvgError(int pageKLLNum, int maxMemoryByte, int SKETCH_BYTE) {
        double dataToSummaryCompression = 1.0*pageN*8/SKETCH_BYTE;
        double SummaryToMemoryCompression = 1.0*pageKLLNum*SKETCH_BYTE/(maxMemoryByte *(1-memoryForMergeBuffer));

        long TOT_SKETCH_N = (long) pageKLLNum *pageN;
        long TOT_SKETCH_SIZE = (long) pageKLLNum*SKETCH_BYTE/8;
        double pageAvgError = 1.0 * TOT_SKETCH_N / TOT_SKETCH_SIZE / 3.0;
        double rate = 1.0 * SKETCH_BYTE * pageKLLNum / (maxMemoryByte *(1-memoryForMergeBuffer));
        long pageStatAvgError;
        if (SummaryToMemoryCompression < 1.0) { // similar to Random Sampling
            pageStatAvgError = (long) Math.ceil(pageAvgError * Math.pow(pageKLLNum, 0.5));
//            if (pageKLLNum <= 10) pageStatAvgError += pageAvgError * 3.0;
        }else /*if(SummaryToMemoryCompression/dataToSummaryCompression<=4.0)*/{
            int memKLLNum = (int)Math.round((maxMemoryByte *(1-memoryForMergeBuffer)) / SKETCH_BYTE);
            long memErr = (long) Math.ceil(pageAvgError * Math.pow(memKLLNum, 0.5));
            pageStatAvgError = (long) Math.ceil(rate * (1.0/3) * memErr + (1-1.0/3) * memErr);
        }/*else
        pageStatAvgError= (long)(approximatePointAvgError(pageKLLNum,(int)(maxMemoryByte *(1-memoryForMergeBuffer)),SKETCH_BYTE)*(1+dataToSummaryCompression/SummaryToMemoryCompression));
*/
        long approxErrPoint = approximatePointAvgError(pageKLLNum,(int)(maxMemoryByte *(1-memoryForMergeBuffer)),SKETCH_BYTE);
        pageStatAvgError=Math.max(pageStatAvgError,approxErrPoint);
        if(SummaryToMemoryCompression/dataToSummaryCompression>=1.0)
            pageStatAvgError = Math.min(pageStatAvgError,(long)(approxErrPoint*(1+dataToSummaryCompression/SummaryToMemoryCompression)));
        return pageStatAvgError;
    }

    private long approximateMergeAvgError(int pageKLLNum, int maxMemoryByte, int SKETCH_BYTE) {
        double dataToSummaryCompression = 1.0*pageN*8/SKETCH_BYTE;
        double SummaryToMemoryCompression = 1.0*pageKLLNum*SKETCH_BYTE/maxMemoryByte;


        long TOT_SKETCH_N = (long) pageKLLNum *pageN;
        long TOT_SKETCH_SIZE = (long) pageKLLNum*SKETCH_BYTE/8;
        double pageAvgError = 1.0 * TOT_SKETCH_N / TOT_SKETCH_SIZE / 3.0;
        double rate = 1.0 * SKETCH_BYTE * pageKLLNum / (maxMemoryByte-SKETCH_BYTE );
        long pageStatAvgError;
        if (SummaryToMemoryCompression < 1.0) { // similar to Random Sampling
            pageStatAvgError = (long) Math.ceil(pageAvgError * Math.pow(pageKLLNum, 0.5));
//            if (pageKLLNum <= 10) pageStatAvgError += pageAvgError * 3.0;
        }else /*if(SummaryToMemoryCompression/dataToSummaryCompression<=4.0)*/{
            int memKLLNum = (maxMemoryByte-SKETCH_BYTE) / SKETCH_BYTE;
            long memErr = (long) Math.ceil(pageAvgError * Math.pow(memKLLNum, 0.5));
            pageStatAvgError = (long) Math.ceil(rate * (1.0/3) * memErr + (1-1.0/3) * memErr);
        }/*else
        pageStatAvgError= (long)(approximatePointAvgError(pageKLLNum,maxMemoryByte-SKETCH_BYTE,SKETCH_BYTE)*(1+dataToSummaryCompression/SummaryToMemoryCompression));
*/
        long approxErrPoint = approximatePointAvgError(pageKLLNum,maxMemoryByte-SKETCH_BYTE,SKETCH_BYTE);
        pageStatAvgError=Math.max(pageStatAvgError,approxErrPoint);
        if(SummaryToMemoryCompression/dataToSummaryCompression>=1.0)
            pageStatAvgError = Math.min(pageStatAvgError,(long)(approxErrPoint*(1+dataToSummaryCompression/SummaryToMemoryCompression)));
        return pageStatAvgError;
    }

    public void testMergeError(int queryPageNum, int maxMemoryByte, int maxSeriesByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time=0, merge_time=0, mergeBuf_time = 0;
        double err_full=0,err_merge=0,err_mergeBuf=0;
        double estimate_err_full=0,estimate_err_merge=0,estimate_err_mergeBuf=0;
        int queryN = queryPageNum*pageN;
        double[] query_a = new double[queryN];
        DoubleArrayList err_record_mergeBuf = new DoubleArrayList();
        DoubleArrayList err_record_merge = new DoubleArrayList();
        DoubleArrayList err_record_full = new DoubleArrayList();

        estimate_err_mergeBuf = 1.0*approximateMergeBufAvgError(queryPageNum,maxMemoryByte,maxSeriesByte)/queryN;
        estimate_err_merge = 1.0*approximateMergeAvgError(queryPageNum,maxMemoryByte,maxSeriesByte)/queryN;
        estimate_err_full = 1.0*approximatePointAvgError(queryPageNum,maxMemoryByte,maxSeriesByte)/queryN;

        for(int T=0;T<TEST_CASE;T++){
            int pageL = random.nextInt(pageNum-queryPageNum+1), pageR = pageL+queryPageNum;
            int posL = pageL*pageN, posR = pageR*pageN;
            double result_full,result_merge;

            merge_time-=new Date().getTime();
            int buf_kll_num = 1;
            List<KLLSketchForQuantile> buf_kll_list = new ArrayList<>(buf_kll_num);
            HeapLongKLLSketch merge_worker = new HeapLongKLLSketch(maxMemoryByte-1*maxSeriesByte);
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


            mergeBuf_time-=new Date().getTime();
            buf_kll_num = (int)(maxMemoryByte/maxSeriesByte*memoryForMergeBuffer);
            buf_kll_list = new ArrayList<>(buf_kll_num);
            HeapLongKLLSketch mergeBuf_worker = new HeapLongKLLSketch((int)(maxMemoryByte*(1-memoryForMergeBuffer)));
            for(int i=pageL;i<pageR;i++){
                if(buf_kll_list.size()<buf_kll_num)
                    buf_kll_list.add(KLLArr[i]);
                else{
//                    System.out.println("\t\t\t\t??!!!!!!!!????\t\t");
//                    mergeBuf_worker.show();
                    mergeBuf_worker.mergeWithTempSpace(buf_kll_list);
//                    mergeBuf_worker.show();
                    buf_kll_list.clear();
                    buf_kll_list.add(KLLArr[i]);
//                    System.out.println("\t\t\t\t??!!!!!!!!????\t\t"+Math.abs(mergeBuf_worker.getN()/2.0-mergeBuf_worker.getApproxRank(doubleToLong(mergeBuf_worker.getN()/2.0))));
                }
            }
            if(!buf_kll_list.isEmpty())
                mergeBuf_worker.mergeWithTempSpace(buf_kll_list);
            mergeBuf_time+=new Date().getTime();


            full_time-=new Date().getTime();
            HeapLongKLLSketch full_worker = new HeapLongKLLSketch(maxMemoryByte);
            for(int i=posL;i<posR;i++)
                full_worker.update(dataToLong(a[i]));
            full_time+=new Date().getTime();

//            if(T==0) {
//                full_worker.show();
//                merge_worker.show();
//                mergeBuf_worker.show();
//                full_worker.showLevelMaxSize();
//                merge_worker.showLevelMaxSize();
//            }

            if (posR - posL >= 0) System.arraycopy(a, posL, query_a, 0, posR - posL);
            Arrays.sort(query_a);

            double q_start=0.01,q_end=0.99,q_add=0.01,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
            for(double q=q_start;q<q_end+1e-10;q+=q_add){
                int query_rank = (int)(q*queryN);

                double mergeBuf_v = longToResult(mergeBuf_worker.findMinValueWithRank(query_rank));
                int mergeBuf_delta_rank = getValueActualRank(query_a,queryN,mergeBuf_v)-query_rank;
                double mergeBuf_relative_err = 1.0*mergeBuf_delta_rank/(queryPageNum*pageN);
                err_mergeBuf+=Math.abs(mergeBuf_relative_err)/(q_count*TEST_CASE);
                err_record_mergeBuf.add(mergeBuf_relative_err);

                double merge_v = longToResult(merge_worker.findMinValueWithRank(query_rank));
                int merge_delta_rank = getValueActualRank(query_a,queryN,merge_v)-query_rank;
                double merge_relative_err = 1.0*merge_delta_rank/(queryPageNum*pageN);
                err_merge+=Math.abs(merge_relative_err)/(q_count*TEST_CASE);
                err_record_merge.add(merge_relative_err);

                double full_v = longToResult(full_worker.findMinValueWithRank(query_rank));
                int full_delta_rank = getValueActualRank(query_a,queryN,full_v)-query_rank;
                double full_relative_err = 1.0*full_delta_rank/(queryPageNum*pageN);
                err_full+=Math.abs(full_relative_err)/(q_count*TEST_CASE);
                err_record_full.add(full_relative_err);

//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
            }
//            double q_start=0.05,q_end=0.95,q_add=0.01,q_count = (q_end-q_start)/q_add+1;
//            for(double q=q_start;q<q_end-1e-10;q+=q_add){
//                int exact_rank = (int)(q*queryN);
//                double v = query_a[exact_rank];
//                long full_rank = full_worker.getApproxRank(v);
//                long merge_rank = merge_worker.getApproxRank(v);
//                err_full+=Math.abs(full_rank-exact_rank)/q_count/TEST_CASE;
//                err_merge+=Math.abs(merge_rank-exact_rank)/q_count/TEST_CASE;
////                System.out.println("?\t\t"+Math.abs(full_rank-exact_rank)+"\t"+Math.abs(merge_rank-exact_rank));
//            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.print("\t\t"+queryPageNum+"\t"+err_mergeBuf+"\t"+err_merge+"\t"+err_full);
        System.out.print("\t\t\tapprox"+"\t"+estimate_err_mergeBuf+"\t"+estimate_err_merge+"\t"+estimate_err_full);
        System.out.println("\t\t\tapprox/exact"+"\t"+estimate_err_mergeBuf/err_mergeBuf+"\t"+estimate_err_merge/err_merge+"\t"+estimate_err_full/err_full);
//        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_mergeBuf-err_full)+"\t"+queryN*(err_merge-err_full));

//        System.out.println("\t\t\t\t\t"+err_merge*(queryPageNum*pageN)+"\t"+err_full*(queryPageNum*pageN));
        String times=fnum.format(1.0*mergeBuf_time/TEST_CASE)+"\t"+fnum.format(1.0*merge_time/TEST_CASE)+"\t"+fnum.format(1.0*full_time/TEST_CASE);
        time_result+=("\t\t"+queryPageNum+"\t"+times)+"\n";
        if(errRecord) {
            err_record_merge.sort(Double::compare);
            err_record_mergeBuf.sort(Double::compare);
            err_record_full.sort(Double::compare);
            File f1 = new File("./" + queryPageNum + "_" + maxMemoryByte + "_" + maxSeriesByte + "_"+dataType+"_KLL");
            FileOutputStream fos1 = new FileOutputStream(f1);
            OutputStreamWriter w1 = new OutputStreamWriter(fos1, StandardCharsets.UTF_8);
            for (int i = 0; i < err_record_merge.size(); i++)
                w1.append("\t").append(String.valueOf(err_record_mergeBuf.getDouble(i))).append("\t").append(String.valueOf(err_record_merge.getDouble(i))).append("\t").append(String.valueOf(err_record_full.getDouble(i))).append("\n");
            w1.close();
            fos1.close();
        }
    }
    public void show_time_result(){
        System.out.println(time_result);
    }
    public static void setTestCase(int tc){TEST_CASE=tc;}




    public static void main(String[] args) throws IOException{
        MainForMergeStatErrorKLL main;
//        main = new MainForMergeStatErrorKLL();
//        main.prepareA();
//        int tmp_seri = 1<<9,tmp_mem=1<<15;
//        main.prepareKLL(tmp_seri);
////        for(int num=1;num<=8;num++)
////            main.testKLL(num,tmp_mem,tmp_seri);
//        for(int num=1;num<=pageNum;num*=2)
//            main.testKLL(num,tmp_mem,tmp_seri);
//        main.show_time_result();

        for(int dataType=0;dataType<=0;dataType++) {
            System.out.println("\n\n\t\tKLL\tdataType:"+dataType+"\t");
            main = new MainForMergeStatErrorKLL();
            MainForMergeStatErrorKLL.dataType = dataType;
            main.prepareA();
            int tmp_seri = 1 << 10, tmp_mem = 1 << 18;
            System.out.println("\tpageN:"+pageN+"\t|summary|:"+tmp_seri+"\t|Memory|:"+tmp_mem);
            main.prepareWorker(tmp_seri);
            for (int num = 1; num <= 256; num *= 2)
                main.testMergeError(num, tmp_mem, tmp_seri);
            for (int num = 512; num <= pageNum; num *= 2) {
//                main.testMergeError(num*15/16, tmp_mem, tmp_seri);
                main.testMergeError(num, tmp_mem, tmp_seri);
            }
            main.show_time_result();
        }
    }
}
