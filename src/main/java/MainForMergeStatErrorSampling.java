import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.eclipse.collections.impl.list.mutable.primitive.LongArrayList;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.util.*;

public class MainForMergeStatErrorSampling {
    static boolean errRecord = true;
    static int dataType = 0;
    static int pageN = 1<<10, pageNum=32768+5000;
    public static int TEST_CASE=16;
    static int N = pageN*pageNum;
    static double[] a,pageMin,pageMax;
    static SamplingHeapForStatMerge[] workerArr;
    String time_result="";

    Random random = new Random(233);

    public void prepareA(){
        a = new double[N];
//        System.out.println("\t\t!"+dataType);
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
        workerArr = new SamplingHeapForStatMerge[pageNum];
        int enoughMemByte = pageN*10;
        for(int i=0;i<pageNum;i++) {
            SamplingHeapForStatMerge worker = new SamplingHeapForStatMerge(enoughMemByte, maxSeriesByte);
            for (int j = 0; j<pageN; j++) worker.update(dataToLong(a[i*pageN+j]));
            worker.compactBeforeSerialization();
            workerArr[i] = worker;
            if(i==0){worker.show();System.out.println();}
        }
    }
    public void preparePageMinMax(){
        pageMin = new double[pageNum];
        pageMax = new double[pageNum];
        for(int i=0;i<pageNum;i++){
            pageMin[i] = Long.MAX_VALUE;
            pageMax[i] = Long.MIN_VALUE;
            double tmp;
            for (int j = 0; j<pageN; j++) {
                tmp=a[i * pageN + j];
                pageMin[i] = Math.min(pageMin[i],tmp);
                pageMax[i] = Math.max(pageMax[i],tmp);
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
    private long dataToLong(double data)  {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }

    private double longToResult(long result)  {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }
    public int calcExactResult(SamplingHeapForStatMerge lastPassWorker,double q,double[] query_a,int queryN,int queryPageNum,double[] queryPageMin,double[] queryPageMax){
        final int K1 = (int)Math.floor(q*queryN), K2 = (int)Math.ceil(q*queryN);
        double cntQ = q;
        LongArrayList lowerBound = lastPassWorker.getLowerBound(cntQ),upperBound = lastPassWorker.getUpperBound(cntQ);
        int cntBoundIndex=0;
        long cntLB =lowerBound.get(cntBoundIndex),cntRB=upperBound.get(cntBoundIndex);

        SamplingHeapForStatMerge worker = new SamplingHeapForStatMerge(lastPassWorker.maxMemByte);
        for(int pass=2;;pass++){
//            if(pass>=8){
//                while(true);
//            }
            System.out.println("\tcntL,R:"+longToResult(cntLB)+"..."+longToResult(cntRB)+"\t\t");
//            System.out.println("\t\t\tcntL,R:"+(cntLB)+"..."+(cntRB)+"\t\tK:"+K1+"  "+K2+"");
//            int FK=0;
//            for(int i=0;i<queryN;i++)if(query_a[i]>=longToResult(cntLB)&&query_a[i]<=longToResult(cntRB))FK++;
//            System.out.println("FK\t\t\t"+FK);
//            FK=0;
//            for(int i=0;i<queryN;i++)if(query_a[i]<longToResult(cntLB))FK++;
//            System.out.println("FK\t\t\t"+FK);
//            FK=0;
//            for(int i=0;i<queryN;i++)if(query_a[i]>longToResult(cntRB))FK++;
//            System.out.println("FK\t\t\t"+FK);

            worker.reset();
            int deltaK=0; // number of data < cntLB.
            long tmp;
            for(int i=0;i<queryPageNum;i++) {
                if(queryPageMin[i]>longToResult(cntRB))continue;
                if (queryPageMax[i] < longToResult(cntLB)) {
                    deltaK+=pageN;
                    continue;
                }
                for (int j = i * pageN, R = j + pageN; j < R; j++) {
                    tmp = dataToLong(query_a[j]);
                    if (tmp >= cntLB && tmp <= cntRB)
                        worker.update(tmp);
                    else if (tmp < cntLB) {
                        deltaK++;
                    }
                }
            }

            if(K1-deltaK<=0||K2-deltaK>worker.totN){
                System.out.println("\t\t\t\t!!!ROLLBACK---------"+"\t\tK1:"+K1+"\tworker_N:"+worker.totN+"\tdeltaK:"+deltaK);
                if(K1-deltaK<=0){
                    cntLB = lowerBound.get(cntBoundIndex+1);
                    cntRB = (K2-deltaK)<=0?(lowerBound.get(cntBoundIndex)-1):(worker.quantile(0.0));
                    cntBoundIndex++;
                    continue;
                }
                else{
                    cntLB = (K1-deltaK>worker.totN)?(upperBound.get(cntBoundIndex) + 1):(worker.quantile(1.0));
                    cntRB = upperBound.get(cntBoundIndex+1);
                    cntBoundIndex++;
                    continue;
                }
            }
            cntQ = 1.0*(K1-deltaK)/ worker.totN;
            System.out.println("\t\t\tafter iteration.\t\tremainK:"+(K1-deltaK)+"  "+(K2-deltaK)+"\tworkerN:"+worker.totN+"\tq:"+cntQ+"sample:"+worker.getSampleSize()+"\tcompression:"+1.0*worker.totN/worker.getSampleSize());
            if(worker.exactResult()){
                long mmp = worker.getSampleK(K1-deltaK-1);

//                System.out.print("\t\tcalc over. pass="+pass);
//                System.out.print("\t\t\t\tresultValue:"+longToResult(mmp));
//                int rk=0;
//                for(int i=0;i<queryN;i++)rk+=query_a[i]<=longToResult(mmp)?1:0;
//                System.out.print("\t\t\t\t"+1.0*rk/queryN+"\t\t"+q+"\n");

                return pass;
            }
//            worker.show();
//            worker.sortSample();for(int i=0;i<worker.indexList.size();i++)System.out.print("\t("+longToResult(worker.value[worker.indexList.get(i)])+")");System.out.println();

            lowerBound = worker.getLowerBound(cntQ);
            upperBound = worker.getUpperBound(cntQ);
            cntBoundIndex = 0;
            cntLB = lowerBound.get(cntBoundIndex);
            cntRB = upperBound.get(cntBoundIndex);
//            System.out.println("\t\t\t\t["+longToResult(lowerBound.get(0))+","+longToResult(upperBound.get(0))
//                +"]\t["+longToResult(lowerBound.get(1))+","+longToResult(upperBound.get(1))+"]\t");
//            int rk=0;
//            for(int i=0;i<queryN;i++)rk+=query_a[i]<longToResult(cntLB)?1:0;
//            System.out.println("\t\t\t\t"+1.0*rk/queryN);
        }

    }


    public void testMergeError(int queryPageNum, int maxMemoryByte, int maxSeriesByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long point_time=0,summary_time=0,full_point_time=0, full_summary_time=0,tmpT=0;
        double err_full=0,err_merge=0;
        int queryN = queryPageNum*pageN;
        double[] query_a = new double[queryN];
        double[] queryPageMin=new double[queryPageNum],queryPageMax=new double[queryPageNum];
        DoubleArrayList err_record_merge = new DoubleArrayList();
        DoubleArrayList err_record_full = new DoubleArrayList();
        double summary_pass = 0.0,point_pass=0.0;

        for(int T=0;T<TEST_CASE;T++){
            int pageL = random.nextInt(pageNum-queryPageNum+1), pageR = pageL+queryPageNum;
            int posL = pageL*pageN, posR = pageR*pageN;
            double result_full,result_merge;
            for(int i=pageL;i<pageR;i++){
                queryPageMin[i-pageL] = pageMin[i];
                queryPageMax[i-pageL] = pageMax[i];
            }


            tmpT = -new Date().getTime();
            SamplingHeapForStatMerge merge_worker = new SamplingHeapForStatMerge(maxMemoryByte);
            for(int i=pageL;i<pageR;i++){
                merge_worker.merge(workerArr[i]);
            }
            tmpT+=new Date().getTime();
            summary_time+=tmpT;
            full_summary_time+=tmpT;


            tmpT = -new Date().getTime();
            SamplingHeapForStatMerge full_worker = new SamplingHeapForStatMerge(maxMemoryByte);
            for(int i=posL;i<posR;i++)
                full_worker.update(dataToLong(a[i]));
            tmpT+=new Date().getTime();
            point_time+=tmpT;
            full_point_time+=tmpT;

//            System.out.println("\t\t\t\t??!!!!!!!!????\t\t"+Math.abs(merge_worker.getN()/2.0-merge_worker.getApproxRank(doubleToLong(merge_worker.getN()/2.0))));
            if(T==0) {
//                merge_worker.show();
//                full_worker.show();
            }
            if (posR - posL >= 0) System.arraycopy(a, posL, query_a, 0, posR - posL);
            Arrays.sort(query_a);
            double q_start=0.1,q_end=0.9,q_add=0.1,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
            for(double q=q_start;q<q_end+1e-10;q+=q_add){

                int query_rank = (int)(q*queryN);
                double merge_v = longToResult(merge_worker.quantile(q));
                int merge_delta_rank = getValueActualRank(query_a,queryN,merge_v)-query_rank;
                double merge_relative_err = 1.0*merge_delta_rank/(queryPageNum*pageN);
                err_merge+=Math.abs(merge_relative_err)/(q_count*TEST_CASE);
                err_record_merge.add(merge_relative_err);

                tmpT = -new Date().getTime();
                summary_pass+=1.0*calcExactResult(merge_worker,q,query_a,queryN,queryPageNum,queryPageMin,queryPageMax)/q_count;
                tmpT+=new Date().getTime();
                full_summary_time+=1.0*tmpT/q_count;


                double full_v = longToResult(full_worker.quantile(q));
                int full_delta_rank = getValueActualRank(query_a,queryN,full_v)-query_rank;
                double full_relative_err = 1.0*full_delta_rank/(queryPageNum*pageN);
                err_full+=Math.abs(full_relative_err)/(q_count*TEST_CASE);
                err_record_full.add(full_relative_err);

                tmpT = -new Date().getTime();
                point_pass+=1.0*calcExactResult(full_worker,q,query_a,queryN,queryPageNum,queryPageMin,queryPageMax)/q_count;
                tmpT+=new Date().getTime();
                full_point_time+=1.0*tmpT/q_count;


//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
            }
        }
        System.out.println("\t\t"+queryPageNum+"\t"+err_merge+"\t"+err_full);
        String times="\tone_pass_time:\t"+fnum.format(1.0*summary_time/TEST_CASE)+"\t"+fnum.format(1.0*point_time/TEST_CASE);
        times+="\tfull_time:\t"+fnum.format(1.0*full_summary_time/TEST_CASE)+"\t"+fnum.format(1.0*full_point_time/TEST_CASE);
        times+="\tpass:\t"+fnum.format(summary_pass/TEST_CASE)+"\t"+fnum.format(point_pass/TEST_CASE);
        time_result+=("\t\t"+queryPageNum+"\t"+times)+"\n";
        if(errRecord) {
            err_record_merge.sort(Double::compare);
            err_record_full.sort(Double::compare);
            File f1 = new File("./" + queryPageNum + "_" + maxMemoryByte + "_" + maxSeriesByte + "_"+dataType+"_Sampling");
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


    public void checkRandom(){
        long START_TIME;

        START_TIME= new Date().getTime();
        XoRoShiRo128PlusRandom random = new XoRoShiRo128PlusRandom(233);
        double sum = 0;int tmp = 100000000;
        for(int i=0;i<tmp;i++)sum+=random.nextDoubleFast();
        System.out.println("\t avg:"+sum/tmp);
        System.out.println("\t\ttime:"+(new Date().getTime()-START_TIME));

        START_TIME= new Date().getTime();
        Random normalRandom = new Random(233);
        sum = 0;
        for(int i=0;i<tmp;i++)sum+=normalRandom.nextDouble();
        System.out.println("\t avg:"+sum/tmp);
        System.out.println("\t\ttime:"+(new Date().getTime()-START_TIME));
    }
    public void checkSampling(){
        long START_TIME = new Date().getTime();
        int maxMemoryByte = 1<<20,maxSeriByte = 1<<15;
        for(int T=1;T<=10;T++) {
            SamplingHeapForStatMerge worker = new SamplingHeapForStatMerge(maxMemoryByte);
            int N = pageN * 1000*T;
            for (int i = 0; i < N; i++)
                worker.update(i);
            worker.sortSample();
            worker.showNum();
            System.out.println("\t\ttime:" + (new Date().getTime() - START_TIME));
            double q_start = 0.01, q_end = 0.99, q_add = 0.01, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
            double avg_err = 0;
            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
                long res = worker.quantile(q);
//            System.out.println("\t\tq:" + q + "\t" +res+"\terr:"+(q*N-res)/N);
                double tmp_err = (q * N - res) / N;
                avg_err += Math.abs(tmp_err) / q_count;
            }
            System.out.println("\t\t\t" + avg_err + "\t\t" + (int) (avg_err * N));
        }
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


//        START_TIME = new Date().getTime();
//        TDigestForStatMerge worker = new TDigestForStatMerge(maxMemoryByte,maxSeriByte);
//        for(int i=0;i<(1<<24);i++)
//            worker.update(i);
//        System.out.println("\t\ttime:"+(new Date().getTime()-START_TIME));
//        START_TIME = new Date().getTime();
//        worker = new TDigestForStatMerge(maxMemoryByte,maxSeriByte);
//        for(int i=0;i<(1<<24);i++)
//            worker.update(i);
//        System.out.println("\t\ttime:"+(new Date().getTime()-START_TIME));


    }

    public static void main(String[] args) throws IOException{
        MainForMergeStatErrorSampling main;

        for(int dataType=0;dataType<=0;dataType++) {
            System.out.println("\n\n\t\tSampling\tdataType:"+dataType+"\t");
            main = new MainForMergeStatErrorSampling();
            MainForMergeStatErrorSampling.dataType = dataType;
            main.prepareA();
            main.preparePageMinMax();
            int tmp_seri = 1 << 9, tmp_mem = 1 << 17;
            main.prepareWorker(tmp_seri);
            for (int num = 128; num <= 16384; num *=2) {
                main.testMergeError(num*3/4, tmp_mem, tmp_seri);
                main.testMergeError(num, tmp_mem, tmp_seri);
            }
            main.show_time_result();
        }
//        main = new MainForMergeStatErrorSampling();
//        main.prepareA();
//        int tmp_seri = 1<<9,tmp_mem=1<<15;
//        main.prepareWorker(tmp_seri);
//        for(int num=1024;num<=pageNum;num*=2)
//            main.testMergeError(num,tmp_mem,tmp_seri);
//        main.show_time_result();
//        main.checkRandom();
//        main.checkSampling();
//        main.checkTDigestSpeed();
    }
}
