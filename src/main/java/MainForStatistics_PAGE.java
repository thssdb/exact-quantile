import org.apache.datasketches.kll.KllDoublesSketch;
import org.apache.iotdb.tsfile.utils.ReadWriteIOUtils;
import org.xerial.snappy.Snappy;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static java.lang.Math.max;

public class MainForStatistics_PAGE {
    static int N = 7000;
    static long[] data = new long[N];
    static long[] data_ordered = new long[N];
    static final long deltaForUnsignedCompare = 1L<<63;

    long seed;
    Random random = new Random(233);
    private long nextLong(){
        return random.nextLong();//&(0xFFFFF0000000000L);
//        seed ^= (seed << 21);
//        seed ^= (seed >>> 35);
//        seed ^= (seed << 4);
//        return seed;
    }
    private long nextLong(int ID,int data_pos){
        long num=-1L;
        switch (ID){
            case 0:
                num = Math.round(1e7 * random.nextGaussian());
                break;
            case 1:
                num = (long) (Math.pow(10, 1 + random.nextDouble() * 17.5)); // log-uniform
                break;
            case 2:
                num = (long) (Math.pow(10, 1 + (17.5 * data_pos / N)));
                break;
            case 3:
                num = (long) (-1e15 + 2e15 * random.nextDouble());
                break;
            case 4:
                num = (long) (-1e15 + 2e15 * data_pos / N);
                break;
            case 5:
                num = (long) (-1e15 + 2e15 * Math.sin(Math.PI * 2 / 1048576 * (data_pos & 1048575)));
                break;
        }
        return num;
    }
    private void reset(){seed=2333L;random.setSeed(233);}
    private void reset(long rnd){seed=rnd;random.setSeed(rnd);}
    private long nextNumber(long x, int i) {
        long num=x;
        num=x&262143;
        if((i&15)<=3)
            num = nextLong();
        return num;
    }
    public void prepare_data(){
//        reset();
        for(int i=0;i<N;i++)data_ordered[i] = data[i] = nextLong();
        Arrays.sort(data_ordered, 0, N);
    }
    public void prepare_data(int DATA_ID){
//        reset();
        for(int i=0;i<N;i++)data_ordered[i] = data[i] = nextLong(DATA_ID, i);
        Arrays.sort(data_ordered, 0, N);
    }
    public double getDelta(long L, long R){
        int K1 = (N+1)/2-1, K2 = N/2;double mid = (data_ordered[K1]+data_ordered[K2])/2.0;
        if(mid<L || mid>R)return 1.0;
        int cnt = 0;
        for(int i=0;i<N;i++) if (data[i]>=L && data[i]<=R) cnt+=1;
        return 1.0*cnt/N;
    }
    public double getDelta(double L, double R){
        int K1 = (N+1)/2-1, K2 = N/2;double mid = (data_ordered[K1]+data_ordered[K2])/2.0;
        if(mid<L || mid>R)return 1.0;
        int cnt = 0;
        for(int i=0;i<N;i++) if (data[i]>=L && data[i]<=R) cnt+=1;
        return 1.0*cnt/N;
    }
    public double getSingleDelta(long X){
        int K1 = (N+1)/2-1, K2 = N/2;
        double mid = (data_ordered[K1]+data_ordered[K2])/2.0;
//        System.out.println("\t\t\t X:"+X+"\t\t\tmid:"+mid);
        int cnt = 0;
        for(int i=0;i<N;i++) if ((data[i]>=X && data[i]<=mid)||(data[i]>=mid && data[i]<=X)) cnt+=1;
        return 1.0*cnt/N;
    }
    public double getSingleDelta(double X){
        int K1 = (N+1)/2-1, K2 = N/2;
        long mid = data_ordered[K1]+(data_ordered[K2]-data_ordered[K1])/2;
//        System.out.println("\t\t\t X:"+X+"\t\t\tmid:"+mid);
        int cnt = 0;
        for(int i=0;i<N;i++) if ((data[i]>=X && data[i]<=mid)||(data[i]>=mid && data[i]<=X)) cnt+=1;
        return 1.0*cnt/N;
    }

    public int getKForKLL(int max_byte) throws IOException{

        int fitK = 0, fitL=10, fitR=max_byte/2;
        while(fitL<fitR){
            int mid = (fitL+fitR+1)/2;
            int cnt_byte = KllDoublesSketch.getMaxSerializedSizeBytes(mid, N, false);
//            KllDoublesSketch worker = KllDoublesSketch.newHeapInstance(mid);
//            for(int i=0;i<N;i++)worker.update(data[i]);
//            byte[] bytes = worker.toByteArray();
//            byte[] com_bytes = Snappy.compress(bytes);
//            int cnt_byte = com_bytes.length;
            if(cnt_byte<=max_byte) {
                fitL = mid;
//                System.out.println("\t\t"+mid+"\t\t"+cnt_byte+"  before_com:"+worker.getSerializedSizeBytes());
//                if(fitL==37) {
//                    System.out.println("K:" + mid + "    " + Arrays.toString(com_bytes));
//                    System.out.println(worker.toString(true, true));
//                    KllDoublesSketch de_worker = KllDoublesSketch.heapify(Memory.wrap(Snappy.uncompress(com_bytes)));
//                    System.out.println(de_worker.toString(true, true));
//                }
            }else fitR = mid-1;
        }
        return fitL;
    }
    public int getKForTDIGEST(int max_byte){return max_byte/12;}

    public List<Double> testKLL_PAGE(int max_byte) throws IOException{
        int KLL_K = getKForKLL(max_byte);
//        System.out.println("KLL K:"+KLL_K);
        KllDoublesSketch worker = KllDoublesSketch.newHeapInstance(KLL_K);
        for(int i=0;i<N;i++)worker.update(data[i]);
        int K1 = (N+1)/2-1, K2 = N/2;
//        System.out.println(worker.getQuantileLowerBound(1.0*K1/N));
//        System.out.println(worker.getQuantileUpperBound(1.0*K1/N));
//        System.out.println(worker.getQuantileLowerBound(1.0*K2/N));
//        System.out.println(worker.getQuantileUpperBound(1.0*K2/N));
        double L = worker.getQuantileLowerBound(1.0*K1/N);
        double R = worker.getQuantileUpperBound(1.0*K2/N);
//        System.out.println("\t\t"+L+","+R);
//        System.out.println("KLL delta:"+getDelta(L,R)*100+"%"+"\t\t theory:"+"(2*)"+worker.getNormalizedRankError(false)*2*100+"%");
//        System.out.println("KLL single delta:"+getSingleDelta(worker.getQuantile(0.5))*100+"%");
        List<Double> final_result = new ArrayList<>(2);
        final_result.add(getDelta(L,R));
        final_result.add(getSingleDelta(worker.getQuantile(0.5)));
        return final_result;
    }

    public List<Double> testTDIGEST_PAGE(int max_byte) throws IOException{
        int TDIGEST_K = getKForTDIGEST(max_byte);
//        System.out.println("???\t\t"+TDIGEST_K);
        TDigestRadixBetterForMedian worker = new TDigestRadixBetterForMedian(TDIGEST_K, TDIGEST_K*6, 64);
        for(int i=0;i<N;i++)worker.add(data[i]);
        int K1 = (N+1)/2-1, K2 = N/2;
        List<Long> result = worker.findResultRange(K1,K2);
        long L = result.get(0), R =result.get(3);
//        System.out.println("\t\t"+L+","+R);
//        System.out.println("TDIGEST delta:"+getDelta(L,R)*100+"%");
//        System.out.println("TDIGEST single delta:"+getSingleDelta(L+(R-L)/2)*100+"%");
        List<Double> final_result = new ArrayList<>(2);
        final_result.add(getDelta(L,R));
        final_result.add(getSingleDelta(getSingleDelta(L+(R-L)/2)));
        return final_result;
    }

    public void test_PAGE(int max_byte)throws IOException {
        reset();
        double KLLAVG = 0, TAVG = 0, KLLWORST = 0, TWORST = 0;
//        double SKLLAVG=0,STAVG=0,SKLLWORST=0,STWORST=0;
        for(int TOT=0;TOT<5;TOT++) {
            reset(nextLong());
            for (int ID = 0; ID < 6; ID++) {
                prepare_data(ID);
//            for(int i=0;i<5;i++)System.out.println("\t\t"+data[i]);
                List<Double> result;

//            result = testKLL_PAGE(max_byte);
//            KLLAVG = (KLLAVG*ID+result.get(0))/(ID+1);
//            KLLWORST = max(KLLWORST, result.get(0));

                result = testTDIGEST_PAGE(max_byte);
                int mmp = TOT*6+ID;
                TAVG = (TAVG * mmp + result.get(0)) / (mmp + 1);
                TWORST = max(TWORST, result.get(0));
            }
        }
//        System.out.println("KLL\t\tAVG:"+KLLAVG*100+"%"+"\t\t\tWORST:"+KLLWORST*100+"%");
        System.out.println("TDIDEST\tAVG:"+TAVG*100+"%"+"\t\t\tWORST:"+TWORST*100+"%");
    }
    public void testSNAPPY() throws IOException {
        final int n = 55;
        byte[] a = new byte[n];
        for(int i=0;i<n;i++)a[i]=(byte)(i%20);
        System.out.println(Arrays.toString(a));
        byte[] a_com = Snappy.compress(a);
        System.out.println(a_com.length);
        System.out.println(Arrays.toString(Snappy.uncompress(a_com)));

        ByteBuffer sb = ByteBuffer.allocate(600);
        int byteLen=0;
        byteLen+=ReadWriteIOUtils.write(2333,sb);
        byteLen+=ReadWriteIOUtils.write(7000,sb);
        byteLen+=ReadWriteIOUtils.write((byte)8,sb);
        byteLen+=ReadWriteIOUtils.write((short)(55), sb);
        for(int i=0;i<11;i++)byteLen+=ReadWriteIOUtils.write((short)(random.nextInt()&31),sb);
        for(int i=0;i<62;i++)byteLen+=ReadWriteIOUtils.write(random.nextLong(),sb);
        System.out.println(byteLen);
        System.out.println(Arrays.toString(sb.array()));
        a_com = Snappy.compress(sb.array());
        System.out.println(Arrays.toString(a_com));
        System.out.println(a_com.length);
        System.out.println(Arrays.toString(Snappy.uncompress(a_com)));
//        ByteBuffer bf = ByteBuffer.allocate(2333);
//        for(int i=0;i<n;i++)bf.put((byte)(i%20));
//        System.out.println(bf.);
//        byte[] b_com = Snappy.compress(bf.toString());
//        System.out.println(b_com.length);
//        System.out.println(Arrays.toString(Snappy.uncompress(b_com)));
    }


    public static void main(String[] args) throws IOException {
        MainForStatistics_PAGE main;
        main = new MainForStatistics_PAGE();
        main.testSNAPPY();
//        main.test_PAGE(512);System.out.println("\n\n");
//        main.test_PAGE(1024);System.out.println("\n\n");
//        main.test_PAGE(2048);System.out.println("\n\n");
    }
}
