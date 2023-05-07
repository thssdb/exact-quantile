import com.google.common.hash.BloomFilter;
import com.google.common.hash.Funnels;
import com.google.common.math.LongMath;
import com.google.common.primitives.Ints;

import java.io.*;
import java.math.RoundingMode;
import java.nio.ByteBuffer;
import java.util.Date;
import java.util.Random;

public class MainForBloomFilter_GUAVA {
    static int N = 100000000;
    final long KK = 47, num_mask = 65536;
    static int TEST_CASE=1;

    long seed;
    Random random = new Random(233);
    private long nextLong(){
        return random.nextLong();
//        seed ^= (seed << 21);
//        seed ^= (seed >>> 35);
//        seed ^= (seed << 4);
//        return seed;
    }
    private void reset(){seed=2333L;random.setSeed(233);}
    private long nextNumber(long x, int i) {
        long num=x;
        num=x&262143;
        if((i&15)<=3)
            num = nextLong();
        return num;
    }


    public long getBFBits(int n,double fpp){
        return (long)((double)(-n) * Math.log(fpp) / (Math.log(2.0D) * Math.log(2.0D)));
    }
    public static double getFPP(double bitsPerKey){
        return Math.exp(-1*bitsPerKey*Math.pow(Math.log(2.0D),2));
    }
    private static long optimalNumOfBits(long n, double p) {
        return (long) (-n * Math.log(p) / (Math.log(2) * Math.log(2)));
    }
    private static int getBFArrayLength(long bits) {
        return Ints.checkedCast(LongMath.divide(bits, 64, RoundingMode.CEILING));
    }
    private static int calcBFSize(long n,double p){
        return getBFArrayLength(optimalNumOfBits(n,p))*8+6;
    }
    public void testBufferToStream() throws IOException {
        ByteBuffer byteBuffer = ByteBuffer.wrap(new byte[233]);
        for(int i=0;i<10;i++)byteBuffer.put((byte)i);
        byteBuffer.rewind();
        System.out.println("\t\t\t" + byteBuffer.get() + "\t" + byteBuffer.get());
        for(int i=0;i<3;i++) {
            InputStream in = new ByteArrayInputStream(byteBuffer.array(), byteBuffer.position()+i*2, 2);
            DataInputStream din = new DataInputStream(in);
            System.out.println("\t\t\t" + din.readByte() + "\t" + din.readByte());
        }
        System.out.println("\t\t\t" + byteBuffer.get() + "\t" + byteBuffer.get());
    }
    public void testGuava() throws IOException {
        int page_N = 8000;double fpp=getFPP(12);
        for(double tmp=1e-2;tmp>1e-7+1e-18;tmp*=0.1)
            System.out.println("\t\t\tfpp:"+tmp+"\tbf bits/key:"+1.0*getBFBits(page_N,tmp)/page_N);
        System.out.println("\t\t\t\tcnt FPP:"+fpp);
//        BloomFilter<Long> bf = BloomFilter.create(Funnels.longFunnel(),page_N,fpp);
//        for(int i=0;i<page_N;i++)bf.put((long) i);
//        int positive=0;
//        for(int i=0;i<1<<20;i++)positive+=bf.mightContain((long)i)?1:0;
//        System.out.println("\t\t\t\tpositive count:"+positive);
//        File f = new File("try_guava.txt");
//        OutputStream of = new FileOutputStream(f);
//        bf.writeTo(of);
//        of.flush();
//        bf=null;
//        InputStream i_f = new FileInputStream(f);
//        bf = BloomFilter.readFrom(i_f,Funnels.longFunnel());
//        positive=0;
//        for(int i=0;i<1<<20;i++)positive+=bf.mightContain((long)i)?1:0;
//        System.out.println("\t\t\t\tpositive count:"+positive);
//        System.out.println("\t\t\t\tcalcSize:"+calcBFSize(page_N,fpp));
    }
    public void testGuavaTime(int bitsPerKey) throws IOException {
        int page_N = 7989;double fpp=getFPP(bitsPerKey);
        BloomFilter<Long> bf = BloomFilter.create(Funnels.longFunnel(),page_N,fpp);
        for(int i=0;i<page_N;i++)bf.put((long) i);
        File f = new File("try_guava.txt");
        OutputStream of = new FileOutputStream(f);
        bf.writeTo(of);
        of.flush();
        bf=null;
        long START_TIME = new Date().getTime();
        for(int T=0;T<TEST_CASE;T++) {
            InputStream i_f = new FileInputStream(f);
            bf = BloomFilter.readFrom(i_f, Funnels.longFunnel());
            for (long i = 0; i < N; i++) bf.mightContain(i);
        }
        long time = new Date().getTime()-START_TIME;
        System.out.println("bits/key:\t"+bitsPerKey+"\t\t\ttime:\t"+time/TEST_CASE+"\t\t\tFPP:"+fpp);
    }
    public void testGuavaBatch(int bitsPerKey) throws IOException {
        int page_N = 7989;double fpp=getFPP(bitsPerKey);
        BloomFilter<Long> bf = BloomFilter.create(Funnels.longFunnel(),page_N,fpp);
        for(int i=0;i<page_N;i++)bf.put((long) i);
    }





    public static void main(String[] args) throws IOException {
        MainForBloomFilter_GUAVA main;
        main = new MainForBloomFilter_GUAVA();
//        main.testBufferToStream();
        main.testGuava();
//        for(int i=10;i<=32;i+=2)
//            main.testGuavaTime(i);

    }
}

// 	[8478]		[1]		[16986]		[55203]		[0]		[33811]
//	(14.29260147472178,81.0)
//	[7872]		[0]		[17326]		[13445]		[88454]
//	(12.564914867994048,76.0)