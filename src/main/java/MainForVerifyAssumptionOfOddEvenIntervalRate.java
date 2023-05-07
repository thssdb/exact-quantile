import it.unimi.dsi.fastutil.ints.IntObjectPair;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.lang3.tuple.MutableTriple;
import org.apache.commons.lang3.tuple.Triple;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

/** Compare sketchByPage & sketchByChunkDivide.
 * Different ChunkSize.
*/
public class MainForVerifyAssumptionOfOddEvenIntervalRate {
    int dataType=-233;
    static int startType=0,endType=3;
    static int pageN = 8192;
    static int N = pageN * 6712/*1000000*/, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 5; // CHECK IT
    public static int MULTI_QUANTILES = 1; // CHECK IT
    static double[] a;
//    static double[] pageMinV,pageMaxV;
    static int compaction_level;
    ArrayList<String> time_result = new ArrayList<>();
    int RESULT_LINE = 0;
    Random random = new Random(233);

    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] = //longToResult(i&4095);
                    random.nextGaussian();
//                    i;
            return;
        }
        BufferedReader reader = null;
        if (dataType == 1||dataType==5)
            reader = new BufferedReader(new FileReader(new File("1_bitcoin.csv")));
        if (dataType == 2||dataType==6)
            reader = new BufferedReader(new FileReader(new File("2_SpacecraftThruster.txt")));
        if (dataType == 3||dataType==7)
            reader = new BufferedReader(new FileReader(new File("3_taxipredition8M.txt")));
        if (dataType == 4||dataType==8)
            reader = new BufferedReader(new FileReader(new File("4_wh.csv")));
        assert reader != null;
        reader.readLine(); // ignore first line.
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if(dataType>=5)a[cntN-1]*=(1.0+random.nextGaussian()*1e-6);
            if (cntN == N) break;
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


    public void testError(int queryN, int maxMemoryByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
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
        double ALLPrSum=0,ALLPrCount=0,FailCount=0,IterCount=0;

        int pointNum=9,pointDelta=queryN/(pointNum+1);
        long[] pointOdd=new long[pointNum],pointEven=new long[pointNum];

        System.out.println("queryN:\t"+queryN+"\t\tMem:"+maxMemoryByte+"\t\tshowPoint:\t"+pointNum+"\t\tTEST_CASE:\t"+TEST_CASE+"\t\tdatasetID:\t"+dataType);

        for (int T = 0; T < TEST_CASE; T++) {
//            System.out.println("\tTEST_ID:"+T);
//            int L = 0, R = queryN;
            int L = LL[T], R = RR[T];
            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);
            ObjectArrayList<MutableTriple<Double,Integer,Integer>>valPosRk=new ObjectArrayList<>(queryN);
            for (int i = L; i < R; i++) valPosRk.add(MutableTriple.of(a[i],i-L,0));
            valPosRk.sort(Comparator.comparingDouble(Triple::getLeft));
            for (int i = 0;i<queryN; i++) valPosRk.get(i).setRight(i);
            valPosRk.sort(Comparator.comparingDouble(Triple::getMiddle));

//            System.out.println("\t\tvalPosRk:\n"+valPosRk);
//
            KLLSketchLazyForVerifyOddEvenIntervalRate sketch = new KLLSketchLazyForVerifyOddEvenIntervalRate(maxMemoryByte);
            for (int i = 0; i < queryN; i++) sketch.update(valPosRk.get(i).getRight());

            if(T==0) {
                System.out.print("\t\tnum_of_compact:\t\t");
                sketch.showCompact();
            }

            int[] oddPreSum=new int[queryN+1],evenPreSum=new int[queryN+1];
            for(IntObjectPair<LongArrayList> compact: sketch.compactRecord){
                LongArrayList values = compact.right();
                long w=1;//1L<<compact.leftInt();
                for(int i=0;i<values.size();i+=2){
                    oddPreSum[(int)values.getLong(i)]+=w;
                    oddPreSum[(int)values.getLong(i+1)]-=w;
                }
                for(int i=1;i<values.size()-1;i+=2){
                    evenPreSum[(int)values.getLong(i)]+=w;
                    evenPreSum[(int)values.getLong(i+1)]-=w;
                }
                evenPreSum[0]+=w;
                evenPreSum[(int)values.getLong(0)]-=w;
                evenPreSum[(int)values.getLong(values.size()-1)]+=w;
                evenPreSum[queryN]-=w;
            }
            for(int i=1;i<queryN;i++){
                oddPreSum[i]+=oddPreSum[i-1];
                evenPreSum[i]+=evenPreSum[i-1];
            }
            long allOdd=0,allEven=0;
            for(int i=0;i<queryN;i++){
                allOdd+=oddPreSum[i];
                allEven+=evenPreSum[i];
            }
            long pointOddSum=0,pointEvenSum=0;
            for(int i=pointDelta,j=0;j<pointNum;i+=pointDelta,j++) {
//                System.out.println("\t\t(" + oddPreSum[i] + " " + evenPreSum[i] + ")");
                pointOddSum+=oddPreSum[i];
                pointEvenSum+=evenPreSum[i];
                pointOdd[j]+=oddPreSum[i];
                pointEven[j]+=evenPreSum[i];
            }
//            System.out.println("\t\tallOdd:\t"+allOdd+"\t\tallEven:\t"+allEven);
//            System.out.println("\t\tpointOdd:\t"+pointOddSum+"\t\tpointEven:\t"+pointEvenSum);
        }

        System.out.println("\tpointID\tpointOdd\tpointEven\t\trelativeOdd");
        for(int i=0;i<pointNum;i++)
            System.out.println("\t"+i+"\t"+pointOdd[i]+"\t"+pointEven[i]+"\t\t"+1.0*pointOdd[i]/(pointOdd[i]+pointEven[i]));
        System.out.println("\n");
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForVerifyAssumptionOfOddEvenIntervalRate main;

        main = new MainForVerifyAssumptionOfOddEvenIntervalRate();
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            main.RESULT_LINE=0;
            main.prepareA(dataType);
//            for (int queryN : new int[]{10000000,20000000})
//                for (int query_mem : new int[]{1024*64,1024*128}) // dataset 0 1 3 表现正常

//            for (int queryN : new int[]{1000000,2000000})
////                for (int query_mem : new int[]{1024*1024,1024*1024*2})
//                for (int query_mem : new int[]{1024*128,1024*256})
//                    main.testError(queryN, query_mem);

//            for (int queryN : new int[]{20000000,40000000})
//                for (int query_mem : new int[]{1024*1024*2,1024*1024*4})
//                    main.testError(queryN, query_mem);
            for (int queryN : new int[]{10000000,20000000})
//                for (int query_mem : new int[]{1024*128,1024*256})
                for (int query_mem : new int[]{1024*1024,1024*1024*2})
                    main.testError(queryN, query_mem);
            System.out.println("\n-------------------------\n");
        }
//        System.out.println("\nError rate:");
//        for (String s : main.err_result)
//            System.out.println(s);
        System.out.println("\t\t\tALL_TIME:"+(new Date().getTime()-START_T));
    }
}
