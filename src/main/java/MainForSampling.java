import it.unimi.dsi.fastutil.doubles.Double2ReferenceAVLTreeMap;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import org.eclipse.collections.api.tuple.primitive.LongDoublePair;
import org.eclipse.collections.impl.tuple.primitive.PrimitiveTuples;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;

public class MainForSampling {
    static boolean errRecord = true;
    static int dataType = 0;
    public static int TEST_CASE=1024;
    static int N = 0;
    static double[] a;
    static long full_time=0;
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


    public void testFullError(int queryN, int maxMemoryByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        full_time=0;
        double err_full=0;
        double[] query_a = new double[queryN];
        DoubleArrayList err_record_full = new DoubleArrayList();
        int FULL_WORKER_SIZE = 0;
        double[] outOfBound = new double[7];

        for(int T=0;T<TEST_CASE;T++){
            int posL = random.nextInt(N-queryN+1), posR = posL+queryN;


            full_time-=new Date().getTime();
            SamplingForStatMerge full_worker = new SamplingForStatMerge(maxMemoryByte);
            for(int i=posL;i<posR;i++)
                full_worker.update(dataToLong(a[i]));
            full_time+=new Date().getTime();


//            System.out.println("\t\t\t\t??!!!!!!!!????\t\t"+Math.abs(merge_worker.getN()/2.0-merge_worker.getApproxRank(doubleToLong(merge_worker.getN()/2.0))));
            if(T==0) {
//                merge_worker.show();
//                full_worker.show();
            }
            if (posR - posL >= 0) System.arraycopy(a, posL, query_a, 0, posR - posL);
            Arrays.sort(query_a);
            double q_start=0.01,q_end=0.99,q_add=0.0001,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
            for(double q=q_start;q<q_end+1e-10;q+=q_add){

                int query_rank = (int)(q*queryN);

                double full_v = longToResult(full_worker.quantile(q));
                int full_delta_rank = getValueActualRank(query_a,queryN,full_v)-query_rank;
                double full_relative_err = 1.0*full_delta_rank/queryN;
                err_full+=Math.abs(full_relative_err)/(q_count*TEST_CASE);
                err_record_full.add(full_relative_err);

//                System.out.println("?\t\tfull:"+full_delta_rank+"\t"+full_relative_err);
                for(int kAvg =1;kAvg<=5;kAvg++)
                    if(Math.abs(full_relative_err)>full_worker.getAvgErr()*kAvg)
                        outOfBound[kAvg-1]+=1.0/(q_count*TEST_CASE);
            }
            FULL_WORKER_SIZE = full_worker.sampleLimit;
        }
//        System.out.println("\t\t"+maxMemoryByte+"\t"+queryN+"\t"+err_full);
        System.out.print(maxMemoryByte+"\tbound:\t\t");for(int i=0;i<outOfBound.length;i++)System.out.print("\t"+outOfBound[i]);System.out.println();
        String times=fnum.format(1.0*full_time/TEST_CASE);
        time_result+=("\t\t"+queryN+"\t"+times)+"\n";
    }
    public void show_time_result(){
        System.out.println("\tTimeResult:");
        System.out.println(time_result);
    }
    public void checkMemCostTreeMap(){
        int NN = 1000000;
        System.gc();
        System.out.println("\t\t\t\tnow mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        long START_MEM = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
        Double2ReferenceAVLTreeMap<LongDoublePair> sampleTreeMap=new Double2ReferenceAVLTreeMap<>();
        for(int i=0;i<NN;i++)sampleTreeMap.put(0.1*i, PrimitiveTuples.pair(0L,0.0d));
        System.gc();
        System.out.println("\t\ttreemap mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()-START_MEM)/NN+"\tBytes for an element");
        System.out.println("\t"+sampleTreeMap.size());
//        sampleTreeMap=null;
//        System.gc();
//        START_MEM = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
//        System.out.println("\t\t\t\tnow mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void checkMemCostHeap(){
        int NN = 1000000;
        System.gc();
        System.out.println("\t\t\t\tnow mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        long START_MEM = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
        SamplingHeapForStatMerge worker_heap=new SamplingHeapForStatMerge(NN*32);
        for(int i=0;i<NN;i++)worker_heap.update(i);
        System.gc();
        System.out.println("\t\t\t\tnow mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        System.out.println("\t\theap mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()-START_MEM)/NN+"\tBytes for an element");

        worker_heap.showNum();
        worker_heap.sortSample();
        System.gc();
        System.out.println("\t\theap mem after sort:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()-START_MEM)/NN+"\tBytes for an element");
        System.out.println("\t\ttest\t"+(worker_heap.indexHeap==null?"null":"not null"));
        System.out.println("\t\t"+worker_heap.quantile(0.5));
        System.gc();
        System.out.println("\t\theap mem after query:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()-START_MEM)/NN+"\tBytes for an element");

        worker_heap.reset();
        System.gc();
        System.out.println("\t\theap mem after reset:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()-START_MEM)/NN+"\tBytes for an element");
System.out.println("\t\ttest\t"+worker_heap.score.length);
    }
    public void checkSpeed(){
        long heapT=0,treeT=0;
        int NN = 100000000,sampleN = 100000;

        SamplingHeapForStatMerge worker_heap = new SamplingHeapForStatMerge(sampleN * 32);
//        SamplingForStatMerge worker_tree = new SamplingForStatMerge(sampleN * 24);
        for(int T=0;T<40;T++) {
            worker_heap.reset();
            heapT -= new Date().getTime();
            for (int i = 0; i < NN; i++) worker_heap.update(i);
            System.out.println("\t\t" + worker_heap.quantile(0.5));
            heapT += new Date().getTime();
//            worker_tree.reset();
//            treeT -= new Date().getTime();
//            for (int i = 0; i < NN; i++) worker_tree.update(i);
//            System.out.println("\t\t" + worker_tree.quantile(0.5));
//            treeT += new Date().getTime();
        }
        System.out.println("\t\t\t"+heapT+"\t\t"+treeT);
    }


    public static void main(String[] args) throws IOException{
        MainForSampling main;
        main = new MainForSampling();
        main.checkSpeed();
//        for(int queryN=640000;queryN<=640000;queryN+=4000){
//            N = queryN*8;
//            main.prepareA();
////            for(int memByte = 24*1;memByte<=24*20;memByte+=24){
////                main.testFullError(queryN,memByte);
////            }
////            for(int memByte = 24*20;memByte<=24*1000;memByte+=24*20){
////                main.testFullError(queryN,memByte);
////            }
//            main.testFullError(queryN,360);
//        }
//        main.checkMemCostTreeMap();
//        main.checkMemCostTreeMap();
//        main.checkMemCostTreeMap();
//        main.checkMemCostHeap();
    }
}
