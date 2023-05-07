import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class MainForTEST {
    static int N = 20000000,K=1;
    final long KK = 47, num_mask = 65536;
    final static int com = 6708;
    static final int mmp = 2;

    long seed;
    Random random = new Random(233);
    private long nextLong(){
        seed ^= (seed << 21);
        seed ^= (seed >>> 35);
        seed ^= (seed << 4);
        return seed;
    }
    private void reset(){seed=2333L;random.setSeed(0);}
    private long nextNumber(long x, int i) {
        long num=x;
//        num=x&262143;
//        num^=1L<<63;
//        num = (x&127)+128;
//        if(Math.abs(nextLong())%10000<=4500)
//            num = nextLong();
        num = (long)(Math.pow(-1,random.nextInt()&1)*Math.pow((random.nextDouble()*(1L<<14)),11));
//        if(random.nextInt(10000)<=1)
//            num = random.nextInt(2000000)-1000000;
//        else num = ((random.nextInt()&1)==0)?Long.MIN_VALUE:Long.MAX_VALUE;
//        if(num<=(-1L<<50))
//            num = (-1L<<50);
//        if(num>=(1L<<50))
//            num = (1L<<50);
//        while(num==Long.MAX_VALUE||num==Long.MIN_VALUE)
//            num = (long)(Math.pow(-1,nextLong()&1)*Math.pow((random.nextDouble()*(1L<<14)),74.0/14));
//        num = (long)(Math.pow(-1,nextLong()&1)*Math.pow((random.nextDouble()*(1L<<16)),64.0/16));
//        num = (long)(Math.pow(-1,nextLong()&1)*Math.exp(43*Math.pow(random.nextDouble(), 0.15)));
        return num;
    }
    public long[] getData(){
        long[] data = new long[N+1];
        reset();
        Random random = new Random(233);
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=/*Math.abs(nextLong())%23+*/1;
            num = nextNumber(num,i);
            data[i-1]=num;
        }
        return data;
    }
    public List<Long> test1(){
        TDigestRadixBetterForMedian worker = new TDigestRadixBetterForMedian(com,6*com, 64);
        reset();
        Random random = new Random(233);
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=/*Math.abs(nextLong())%23+*/1;
            num = nextNumber(num,i);
            worker.add(num);
        }
        List<Long> result = new ArrayList<>(mmp*2);
        for(int i=0;i<mmp;i+=2)
            result.addAll(worker.findResultRange((long)(1.0*N/mmp*i), (long)(1.0*N/mmp*(i+1))));
//        System.out.println(result);
        return result;
    }
    public List<Long> test2(){
        EclipseCollectionsHashMapForQuantile worker = new EclipseCollectionsHashMapForQuantile(64,16);
        reset();
        Random random = new Random(233);
        for(int i=1;i<=N;i++) {
            long num = nextLong(), freq =/*Math.abs(nextLong())%23+*/1;
            num = nextNumber(num, i);
            worker.insert(num, freq);
        }
        List<Long> result = new ArrayList<>(mmp*2);
        for(int i=0;i<mmp;i+=2) {
            List<Long> temp = worker.findResultIndex((long)(1.0*N/mmp*i), (long)(1.0*N/mmp*(i+1)));
            long bits = worker.bitsOfValue-worker.getRemainingBits();
            result.add(temp.get(0)<<bits);
            result.add((temp.get(0)<<bits)-1+(1L <<bits));
            result.add(temp.get(2)<<bits);
            result.add((temp.get(2)<<bits)-1+(1L <<bits));
        }
        return result;
    }

    public static void signedShowAbsoluteRankError(long[] a, List<Long> result){
        DecimalFormat df1 = new DecimalFormat("#0.0");
        DecimalFormat df2 = new DecimalFormat("#.#E0");
        for(int i=1;i<mmp;i++)
            System.out.print("\t"+df1.format(1.0/mmp*i)+"\t\t");
        System.out.println();
        for(int i=1;i<mmp;i++){
            long L = result.get(i*2);
            long R = result.get(i*2+1);
            long expectedRank = (long)(1.0*N/mmp*i);
            long actualRank = 0;
            long rest = 0;
            for(int j=0;j<N;j++)
                if(a[j]>=L&&a[j]<=R)
                    rest++;
            System.out.print("\t"+/*df2.format(1.0*rest/N)*/rest+"");
            System.out.print(" "+L+" "+R+"("+(R-L)+")\t");
        }
        System.out.println();
    }

    public static void unsignedShowAbsoluteRankError(long[] a, List<Long> result){
        long deltaForUnsigned = 1L<<63;
        int mmp = 10;
        for(int i=1;i<mmp;i++)
            System.out.print("\t"+1.0/mmp*i+"\t");
        System.out.println();
        for(int i=1;i<mmp;i++){
            long L = result.get(i*2)^deltaForUnsigned;
            long R = result.get(i*2+1)^deltaForUnsigned;
            long expectedRank = (long)(1.0*N/mmp*i);
            long actualRank = 0;
            long rest = 0;
            for(int j=0;j<N;j++)
                if((a[j]^deltaForUnsigned)>=L&&(a[j]^deltaForUnsigned)<=R)
                    rest++;
            System.out.print("\t"+(1.0*rest/N)+"\t");
        }
        System.out.println();
    }

    public static void testBufferToStream() throws IOException {
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
    public static void testByteBuffer()throws IOException {
        ByteBuffer byteBuffer = ByteBuffer.wrap(new byte[233]);
        for(int i=0;i<30;i++)byteBuffer.put((byte)i);
        byteBuffer.rewind();
        System.out.println("\t\t\t" + byteBuffer.get() + "\t" + byteBuffer.get());
        ByteBuffer byteBuffer1 = byteBuffer.slice();
    }

    public static void main(String[] args) {

        MainForTEST main = new MainForTEST();
        List<Long> result;
        long[] data = main.getData();

        result = main.test1();
        signedShowAbsoluteRankError(data, result);
//        for(int i=0;i<20;i++)System.out.println(data[i]);
        int numMIN=0, numMAX=0;
        for(int i=0;i<N;i++)
            if(data[i]==Long.MIN_VALUE)numMIN++;
        else if(data[i]==Long.MAX_VALUE)numMAX++;
        System.out.println("\t\t\t numMINMAX:"+numMIN+"  "+numMAX);
    }
}
