import java.util.Random;

public abstract class TDigestForPageStat {
  protected static final int BUCKET_K = 3; // 2^BUCKET_K
  long XORSHIFT=new Random().nextInt();//0x2333333319260817L;
//  Random test_random = new Random();

  protected short getNextRand(int toAnd){ // xor shift *
    XORSHIFT^=XORSHIFT>>>12;
    XORSHIFT^=XORSHIFT<<25;
    XORSHIFT^=XORSHIFT>>>27;
    return (short) ((XORSHIFT*0x2545F4914F6CDD1DL)&toAnd);
  }
  protected int getNextRand(){ // xor shift *
    XORSHIFT^=XORSHIFT>>>12;
    XORSHIFT^=XORSHIFT<<25;
    XORSHIFT^=XORSHIFT>>>27;
    return (int) ((XORSHIFT*0x2545F4914F6CDD1DL));
  }

  public TDigestForPageStat() {
  }

}