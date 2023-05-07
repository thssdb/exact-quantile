/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */

import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.eclipse.collections.api.map.primitive.MutableLongLongMap;
import org.eclipse.collections.api.tuple.primitive.LongLongPair;
import org.eclipse.collections.impl.map.mutable.primitive.LongLongHashMap;

import java.util.ArrayList;
import java.util.List;

/**
 * hppc long-long hashmap with <= 65536 elements.
 * rebuild (right shift keys and merge) to keep <=65536.
 */
public class GCHashMapForMem {
  private final int maxSize;
  int bitsOfValue, remainingBits, minBits;

  MutableLongLongMap hashMap,newMap;
  int totSize;
  long[] value, count;
  private long deltaForUnsignedCompare;
  IntArrayList index;
  long DEBUG;


  public GCHashMapForMem(int bits, int minBits, int maxSize) {
    this.minBits = minBits;
    bitsOfValue = remainingBits = bits;
    if (bits == 64) deltaForUnsignedCompare = 1L << 63; // unsigned long
    else deltaForUnsignedCompare = 0;
    value = new long[maxSize];
    count = new long[maxSize];
    hashMap = new LongLongHashMap(maxSize);
    this.maxSize = maxSize;
    index = new IntArrayList(maxSize);
    newMap = new LongLongHashMap(maxSize);
  }

  private void rebuild() { // called when total size == maxSize
//    System.out.println("[rebuild]+remaining:"+remainingBits+"  tot:"+totSize);
    int SHR = 1;
    deltaForUnsignedCompare = 0;
    newMap.clear();
    hashMap.forEachKeyValue((k,v)->newMap.addToValue(k>>>1,v));

    while (newMap.size() == maxSize) {
      SHR++;
      newMap.clear();
//      for(LongLongPair p:hashMap.keyValuesView())
//        newMap.addToValue(p.getOne()>>>SHR,p.getTwo());
      final int shr = SHR;
      hashMap.forEachKeyValue((k,v)->newMap.addToValue(k>>>shr,v));
    }
    remainingBits -= SHR;
    hashMap = newMap;
  }

  public void insert(long num, long freq) {
    num >>>= bitsOfValue - remainingBits;
    hashMap.addToValue(num,freq);
    if (hashMap.size() == maxSize)
      rebuild();
  }
  public void insertOne(long num) {
    num >>>= bitsOfValue - remainingBits;
    hashMap.addToValue(num,1L);
    if (hashMap.size() == maxSize)
      rebuild();
  }

  public int getRemainingBits() {
    return remainingBits;
  }


  public List<Long> findResultIndex(long K1, long K2) {
    List<Long> result = new ArrayList<>(8);
    long sum = 0;

    totSize = hashMap.size();
    index.size(totSize);
    for (int i = 0; i < totSize; i++) index.set(i, i);
    int tmp = 0;

    for (LongLongPair p : hashMap.keyValuesView()) {
      value[tmp] = p.getOne();
      count[tmp] = p.getTwo();
      tmp++;
    }
    index.sort((x, y) -> Long.compare(value[x] ^ deltaForUnsignedCompare, value[y] ^ deltaForUnsignedCompare));
    int x;
    for (int i = 0; i < totSize; i++) {
      x = index.getInt(i);
//      System.out.println(count[x] + "  " + value[x]);
      sum += count[x];
      if (sum >= K1 && result.size() == 0) {
        result.add(value[x]);
        result.add(sum - count[x]);
      }
      if (sum >= K2 && result.size() == 2) {
        result.add(value[x]);
        result.add(sum - count[x]);
        break;
      }
    }

    return result;
  }




  public void reset(int bits, int minBits){
    remainingBits = bits;
    this.minBits = minBits;
    hashMap.clear();
  }
}
