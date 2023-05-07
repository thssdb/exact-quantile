import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;

/** Compare sketchByPage & sketchByChunkDivide.
 * Different ChunkSize.
*/
public class MainForExactTDigest {
    int dataType=-233;
    static int startType=0,endType=0;
    static int pageN = 8192;
    static int N = pageN * 6712*2/*6712*/, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 1,QUANTILE_PER_TEST=5; // CHECK IT
    static double[] a;
//    static double[] pageMinV,pageMaxV;
    ArrayList<String> time_result = new ArrayList<>();
    int RESULT_LINE = 0;
    Random random = new Random(233);

    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        LogNormalDistribution Rlog1 = new LogNormalDistribution(10,0.001);
        LogNormalDistribution Rlog2 = new LogNormalDistribution(5,0.1);
        UniformRealDistribution R01 = new UniformRealDistribution(0,1);
        double Emax=300;
        //    DoubleArrayList data = new DoubleArrayList();
        //    for(int i=0;i<(int)1e3;i++){
        //        double x = Math.pow(10, Emax*(Math.pow(Rlog1.sample(),2)*2-1));
        //        if(bit.nextBoolean())
        //            x=-x;
        //        data.add(x);
        //    }
        System.out.println("\t\t"+Double.MAX_VALUE);
        if (dataType == 0) {
            for (int i = 0; i < N; i++) {
                a[i] = //longToResult(i&4095);
                    //Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
//                    random.nextGaussian();
                    Rlog1.sample();
//                    i-0.5*Math.floor(i/pageN)*pageN;
//                    i;
                if(R01.sample()<0.1)
                    a[i] = Math.pow(10, Emax*(Math.pow(R01.sample(),2)*2-1));
//                if(R01.sample()<0.5)
//                    a[i]=Rlog2.sample();
            }
            return;
        }
        BufferedReader reader = null;
        if (dataType == 1)
            reader = new BufferedReader(new FileReader(new File("1_bitcoin.csv")));
        if (dataType == 2)
            reader = new BufferedReader(new FileReader(new File("2_SpacecraftThruster.txt")));
        if (dataType == 3)
            reader = new BufferedReader(new FileReader(new File("3_taxipredition8M.txt")));
        if (dataType == 4)
            reader = new BufferedReader(new FileReader(new File("4_wh.csv")));
        assert reader != null;
        reader.readLine(); // ignore first line.
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if (cntN == N) break;
        }
    }

//    public void prepareWorker(int maxSeriesByte) {
//        pageMinV = new double[pageNum];pageMaxV = new double[pageNum];
//
//        int enoughMemByte = pageN * 9;
//        for (int i = 0; i < pageNum; i++) {
//            pageMinV[i] = Long.MAX_VALUE;
//            pageMaxV[i] = Long.MIN_VALUE;
//            for (int j = 0; j < pageN; j++) {
//                double v =a[i*pageN+j];
//                pageMinV[i] = Math.min(pageMinV[i],v);
//                pageMaxV[i]=Math.max(pageMaxV[i],v);
//            }
////            worker.show();
//        }
//    }

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

    public void testError(int queryN, int maxMemoryByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        long full_time = 0, merge_page_time = 0;
        double avg_iteration=0;
        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }

        for (int T = 0; T < TEST_CASE; T++) {
//            System.out.println("\tTEST_ID:"+T);
//            int L = 0, R = queryN;
            int L=LL[T],R=RR[T];

//            int MMP_COUNT=0;
//            for(int i=L;i<R;i++)if(a[i]<-74.00054931641387)MMP_COUNT++;
//            System.out.println("[][][]MMP][][]"+MMP_COUNT);

            merge_page_time -= new Date().getTime();
            TDigestForExact tDigestWorker = new TDigestForExact(maxMemoryByte);


            for (int i = L; i < R; i++)
                tDigestWorker.update(a[i]);
//            tDigestWorker.show();
//System.out.println("?????1");
            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);

            int q_count=QUANTILE_PER_TEST;
            double ratioPerQuery = 1.0 / (q_count * TEST_CASE);
            for(int INNER_T=0;INNER_T<q_count;INNER_T++)
            {
                double q = random.nextDouble();
                int query_rank1 = (int) Math.floor(q * (queryN-1)+1),query_rank2 = (int) Math.ceil(q * (queryN-1)+1);
                long last_n=queryN;

                double[] iterate_result;// = tDigestWorker.findResultRange(query_rank1,query_rank2);
                iterate_result =tDigestWorker.getFilter(0,0,0,0,query_rank1,query_rank2);
                last_n = tDigestWorker.findMaxNumberInRange(iterate_result[0], iterate_result[1]);
//                System.out.println("?????2");
                avg_iteration+=ratioPerQuery;
                int MMP=0;
                while(iterate_result[0]<iterate_result[1]&&iterate_result.length==2) {
                    if(++MMP>20)break;
                    avg_iteration+=ratioPerQuery;
                    TDigestForExact cntWorker = new TDigestForExact(maxMemoryByte);

                    double valL=iterate_result[0],valR=iterate_result[1];
                    int CountOfLessThanValL = 0,CountOfValL=0,CountOfValR=0;
                    System.out.println("\t(start TD Iteration:"+"\t\tvalL,R:"+valL+","+valR+"\t\t\tMMP:\t"+MMP+"\t\tn:"+last_n);
                    if(last_n<=maxMemoryByte/8){
                        System.out.println("\t\t\tlastPass\t\tmemUsedRate:\t"+1.0*last_n*8/maxMemoryByte+"\t\tlast_n:\t"+last_n);
                        DoubleArrayList tmpList = new DoubleArrayList((int)last_n);
                        for(int i=L;i<R;i++) {
                            if (a[i] >= valL && a[i] <= valR) tmpList.add(a[i]);
                            else if (a[i] < valL) CountOfLessThanValL++;
                        }
                        if(tmpList.size()>maxMemoryByte/8)System.out.println("!!!!!!!FAILED ALL VALUE TOO MUCH");
                        int k1=query_rank1-CountOfLessThanValL-1,k2=query_rank2-CountOfLessThanValL-1;
//                        System.out.println("\t\t\t\t use brute force. k1,k2:\t"+k1+","+k2+"\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tquery_rank:"+query_rank1+","+query_rank2);
                        tmpList.sort(Double::compare);
                        iterate_result[0] = tmpList.getDouble(k1);
                        iterate_result[1] = tmpList.getDouble(k2);
                        break;
                    }else {
                        for (int i = L; i < R; i++) {
                            if (a[i] > valL && a[i] < valR) cntWorker.update(a[i]);
                            else if (a[i] < valL) CountOfLessThanValL++;
                            else if(a[i]==valL)CountOfValL++;
                            else if(a[i]==valR)CountOfValR++;
                        }
//                        cntWorker.show();
//                        System.out.println("\t\t\t\t\t\tCountOfLessThanValL:" + CountOfLessThanValL + "\t\tcntSketch_N:" + cntWorker.totN);
                        int cntRank1 = query_rank1 - CountOfLessThanValL;
                        int cntRank2 = query_rank2 - CountOfLessThanValL;
//                        System.out.println("\t\t\t\t\t\tcntRank:" + cntRank1 + " " + cntRank2);
//                        iterate_result = cntWorker.findResultRange(cntRank1, cntRank2);
                        iterate_result =cntWorker.getFilter(CountOfValL,CountOfValR,valL,valR,cntRank1,cntRank2);
                        iterate_result[0] = Math.max(iterate_result[0], valL);
                        iterate_result[1] = Math.min(iterate_result[1], valR);
                        last_n = cntWorker.findMaxNumberInRange(iterate_result[0], iterate_result[1]);
//                        System.out.println("\t\t\t\tcntL,R:" + iterate_result[0] + "," + iterate_result[1] + "\t\t\tCountOfLessThanValL:" + CountOfLessThanValL + "\t\tcntRank:" + cntRank1 + " " + cntRank2+"\t\tlastN:"+last_n);
                    }
                }
                double exact_quantile_v = (iterate_result[0]+iterate_result[1])*0.5;
                System.out.println("\t\tq:"+q+"\texact_quantile="+exact_quantile_v);
            }
        }
        System.out.println("\t\t"+avg_iteration);
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
//        System.out.println("TEST_CASE="+TEST_CASE+"\tDATASET:"+dataType+"\tqueryN:\t" + queryN+"\tmemory:\t" + maxMemoryByte+"\t\tavg_iteration:\t"+avg_iteration);
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForExactTDigest main;

        main = new MainForExactTDigest();
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            System.out.println("\t|||DATASET:"+"\t"+dataType);
            main.RESULT_LINE=0;
            main.prepareA(dataType);
            int[] nn=new int[]{100000,400000,100000000};
            int[][] mm=new int[][]{new int[]{1024*12,1024*24,1024*32},new int[]{1024*32,1024*48,1024*64},new int[]{1024*512}};
            int[] ttcc=new int[]{10,2,1};

            for(int ii=2;ii<3;ii++)
                for(int jj=0;jj<1;jj++){
                    int queryN=nn[ii],queryMem=mm[ii][jj];
                    TEST_CASE*=ttcc[ii];
                    System.out.println("\t\tqueryN:\t"+queryN+"\tMemory\t"+queryMem+"\t");
                    System.out.println("\t\tavgPass");

                    main.testError(queryN, queryMem);
                    TEST_CASE/=ttcc[ii];
                }
        }
//        System.out.println("\nError rate:");
//        for (String s : main.err_result)
//            System.out.println(s);
        System.out.println("\t\t\tALL_TIME:"+(new Date().getTime()-START_T));
    }
}
