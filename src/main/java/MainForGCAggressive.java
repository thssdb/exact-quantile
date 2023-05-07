import java.util.Date;
import java.util.List;
import java.util.Random;

public class MainForGCAggressive {
    int N = 20000000;
    final long KK = 47, num_mask = 65536;
    final static int com = 6708;

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
//        num&=2147483647;
//        if(random.nextInt(10000)<=10)
//            num = random.nextInt(2000000);
//        else num = ((random.nextInt()&1)==0)?0:-1L;
        num = (long)(Math.pow(-1,random.nextInt()&1)*Math.pow((random.nextDouble()*(1L<<14)),11));
//        num = (x&127)+128;
//        if(Math.abs(nextLong())%10000<=6000)
//            num = (long)(random.nextGaussian()*1e17);
        return num;
    }
    public void test1(){
//        System.out.println("   start mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        EclipseCollectionsHashMapForQuantile worker = new EclipseCollectionsHashMapForQuantile(64,16);
        reset();
        Random random = new Random(233);
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=/*Math.abs(nextLong())%23+*/1;
            num = nextNumber(num,i);
            worker.insert(num, freq);
        }
        List<Long> result = worker.findResultIndex(N/2, N/2+1);
        System.out.println("\t\t\t\t"+result.toString());
//        System.out.println(worker.DEBUG);
        System.out.println(worker.getRemainingBits());
//        System.out.println("   end mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }
    public void test2(){
//        System.out.println("   start mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
        GCHashMapComplexForQuantile worker = new GCHashMapComplexForQuantile(64,16, 0.5);
        reset();
        Random random = new Random(233);
        for(int i=1;i<=N;i++) {
            long num=nextLong(),freq=/*Math.abs(nextLong())%23+*/1;
            num = nextNumber(num,i);
            worker.insert(num, freq);
        }
        List<Long> result = worker.findResultIndex(N/2, N/2+1);
        System.out.println("\t\t\t\t"+result.toString());
        System.out.println(worker.DEBUG);
        System.out.println(worker.getRemainingBits());
//        System.out.println("   end mem:"+1.0*(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    }

    public static void main(String[] args) {

        int TEST_CASE = 1, sum;
        long START_TIME = new Date().getTime();


        sum=0;
        for(int i=0;i<TEST_CASE;i++) {
            START_TIME = new Date().getTime();
            MainForGCAggressive main = new MainForGCAggressive();
            main.test1();
            sum+=(new Date().getTime() - START_TIME);
        }
        System.out.println("gc: "+sum/TEST_CASE+"ms");


        sum=0;
        for(int i=0;i<TEST_CASE;i++) {
            START_TIME = new Date().getTime();
            MainForGCAggressive main = new MainForGCAggressive();
            main.test2();
            sum+=(new Date().getTime() - START_TIME);
        }
        System.out.println("gc aggressive: "+sum/TEST_CASE+"ms");
    }
}
