import it.unimi.dsi.fastutil.doubles.DoubleArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class EvaluatingMergingKLLSketch {
    int dataType = 1;
    static int pageN = 1<<12, pageNum=/*32768*2+5000*/13000;
    public static int TEST_CASE=64;
    static int N = pageN*pageNum;
    static double[] a;
    static KLLSketchForQuantile[] KLLArr;
    String time_result="";

    Random random = new Random(233);

    public void prepareA(int dataType) throws IOException {
        if(a==null)a = new double[N];
        this.dataType = dataType;

        if(dataType==0){
            for(int i=0;i<N;i++)a[i]=Math.pow(-1,random.nextInt(2))*Math.pow(10.0,(2*Math.pow(random.nextDouble(),2)-1)*300);
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
            BufferedReader reader = new BufferedReader(new FileReader(new File("2_physiological_stress.txt")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN= 0;
            while((line=reader.readLine())!=null){
                a[cntN++] = Double.parseDouble(line);
                if(cntN==N)return;
            }
        }
        if(dataType==3){
            BufferedReader reader = new BufferedReader(new FileReader(new File("4_taxipredition8M.txt")));
            reader.readLine(); // ignore first line.
            String line;
            int cntN= 0;
            while((line=reader.readLine())!=null){
                a[cntN++] = Double.parseDouble(line);
                if(cntN==N)return;
            }
        }
        if(dataType==4){
            BufferedReader reader = new BufferedReader(new FileReader(new File("5_wh.csv")));
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

    public void testMergeError(int queryPageNum, int maxMemoryByte, int maxSeriesByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time=0, merge_time=0;
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
            int buf_kll_num = 1;
            List<KLLSketchForQuantile> buf_kll_list = new ArrayList<>(buf_kll_num);
            HeapLongStrictKLLSketch merge_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
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



            full_time-=new Date().getTime();
            HeapLongStrictKLLSketch full_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            for(int i=posL;i<posR;i++)
                full_worker.update(dataToLong(a[i]));
            full_time+=new Date().getTime();

            if(T==0) {
//                merge_worker.show();
//                full_worker.show();
//                merge_worker.showLevelMaxSize();
//                full_worker.showLevelMaxSize();
            }

            if (posR - posL >= 0) System.arraycopy(a, posL, query_a, 0, posR - posL);
            Arrays.sort(query_a);

            double q_start=0.01,q_end=0.99,q_add=0.005,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
            for(double q=q_start;q<q_end+1e-10;q+=q_add){
                int query_rank = (int)(q*queryN);

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
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("\t\t"+queryPageNum+"\t"+err_merge+"\t"+err_full);//        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_mergeBuf-err_full)+"\t"+queryN*(err_merge-err_full));

//        System.out.println("\t\t\t\t\t"+err_merge*(queryPageNum*pageN)+"\t"+err_full*(queryPageNum*pageN));
        String times="\t"+fnum.format(1.0*merge_time/TEST_CASE)+"\t"+fnum.format(1.0*full_time/TEST_CASE);
        time_result+=("\t\t"+queryPageNum+"\t"+times)+"\n";
//        if(errRecord) {
//            err_record_merge.sort(Double::compare);
//            err_record_full.sort(Double::compare);
//            File f1 = new File("./" + queryPageNum + "_" + maxMemoryByte + "_" + maxSeriesByte + "_"+dataType+"_KLL");
//            FileOutputStream fos1 = new FileOutputStream(f1);
//            OutputStreamWriter w1 = new OutputStreamWriter(fos1, StandardCharsets.UTF_8);
//            for (int i = 0; i < err_record_merge.size(); i++)
//                w1.append("\t").append(String.valueOf(err_record_merge.getDouble(i))).append("\t").append(String.valueOf(err_record_full.getDouble(i))).append("\n");
//            w1.close();
//            fos1.close();
//        }
    }
    public void show_time_result(){
        System.out.println(time_result);
    }
    public static void setTestCase(int tc){TEST_CASE=tc;}




    public static void main(String[] args) throws IOException{
        EvaluatingMergingKLLSketch main;
//        main = new MainForMergeStatErrorKLL();
//        main.prepareA();
//        int tmp_seri = 1<<9,tmp_mem=1<<15;
//        main.prepareKLL(tmp_seri);
////        for(int num=1;num<=8;num++)
////            main.testKLL(num,tmp_mem,tmp_seri);
//        for(int num=1;num<=pageNum;num*=2)
//            main.testKLL(num,tmp_mem,tmp_seri);
//        main.show_time_result();

        for(int dataType=4;dataType<=4;dataType++) {
            System.out.println("\n\n\t\tKLL\tdataType:"+dataType+"\t");
            main = new EvaluatingMergingKLLSketch();
            main.prepareA(dataType);
            int tmp_seri = 1 << 10, tmp_mem = 1 << 17;
            System.out.println("\tpageN:"+pageN+"\t|summary|:"+tmp_seri+"\t|Memory|:"+tmp_mem);
            main.prepareWorker(tmp_seri);
//            for (int num = 1; num <= 256; num *= 2)
//                main.testMergeError(num, tmp_mem, tmp_seri);
//            for (int num = 32; num <= pageNum; num *= 2) {
//                main.testMergeError(num, tmp_mem, tmp_seri);
//            }
//            main.testMergeError(12207, tmp_mem, tmp_seri); // 5E7
            for (int num = 512; num <= 512+1024*8; num +=1024) {
                main.testMergeError(num, tmp_mem, tmp_seri);
                if(num==512)
                    main.testMergeError(1024, tmp_mem, tmp_seri);
            }
            main.show_time_result();
        }
    }
}
