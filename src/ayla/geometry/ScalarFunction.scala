/**
 * Ayla Geometry Core
 * (c) 2011-2014 William Harvey
 * http://www.aylasoftware.org
 */
package ayla.geometry

import java.util.Comparator
import scala.collection.JavaConversions._
import ayla.util.UnionFind
import scala.annotation.tailrec

class ScalarFunction(val sc: SimplicialComplex, val getFuncVal: Array[Float]) {
  val minFuncVal = sc.vertices.indices.map(getFuncVal).min
  val maxFuncVal = sc.vertices.indices.map(getFuncVal).max
  
  val rangeFuncVal = if (maxFuncVal - minFuncVal != 0) maxFuncVal - minFuncVal else 1

  class VertComparator extends java.util.Comparator[Int] with Serializable {
    override def compare(v1: Int, v2: Int): Int = {
      val f1 = getFuncVal(v1)
      val f2 = getFuncVal(v2)

      if (f1 < f2) {
        return -1
      } else if (f1 > f2) {
        return 1
      } else if (v1 < v2) {
        return -1
      } else if (v1 > v2) {
        return 1
      } else {
        return 0
      }
    }
  }

  val vc = new VertComparator

  def getGradientMag(v1: Int, v2: Int): Double = {
    val d = dist(v1, v2)
    if (d == 0) {
      return 0
    }
    return math.abs(getFuncVal(v1) - getFuncVal(v2)) / d
  }

  def dist(v1: Int, v2: Int): Double = {
    val nDims = sc.vertices(0).length
    
    @tailrec
    def getDist(i: Int = 0, sumSq: Double = 0): Double = {
    	if (i == nDims)
    	  math.sqrt(sumSq)
    	else {
    	  val diff = sc.vertices(v1)(i) - sc.vertices(v2)(i)
    	  getDist(i + 1, sumSq + (diff * diff))
    	}
    }
    
    getDist()
  }

  def getLargerNeighborsOf(v: Int): Iterator[Int] = sc.getNeighbors(v).iterator.filter(n => vc.compare(v, n) < 0)
  def getSmallerNeighborsOf(v: Int): Iterator[Int] = sc.getNeighbors(v).iterator.filter(n => vc.compare(v, n) > 0)

  class StandardPersistencePair(extremum: Int, saddle: Int, killedBy: Int) extends StandardPersistencePairLike(extremum, saddle, killedBy) {
    override val persistence = math.abs(getFuncVal(extremum) - getFuncVal(saddle))
  }

  def getStandardPersistencePairs(): Array[this.StandardPersistencePair] = {
    val n = sc.vertices.size
    val filtration = (0 until n).toArray.sortWith(vc.compare(_, _) < 0)
    val pairsSweepUp = sweep(filtration, getSmallerNeighborsOf, vc.compare(_, _))
    val pairsSweepDown = sweep(filtration.reverse, getLargerNeighborsOf, -vc.compare(_, _))

    return (pairsSweepUp ::: pairsSweepDown).toArray
  }

  private def sweep(vertFiltration: Array[Int], link: Int => Iterator[Int], comp: (Int, Int) => Int): List[StandardPersistencePair] = {
    val n = sc.vertices.size
    val uf = new UnionFind(n) //new HashUnionFind((0 until n).toSet)
    val componentRep = new scala.collection.mutable.HashMap[Int, Int]

    val mostExtreme = vertFiltration.reduceLeft((x, y) => if (comp(x, y) < 0) x else y)
    println("Global:  " + mostExtreme)

    val extrema = new scala.collection.mutable.HashSet[Int]
    val pairs = new scala.collection.mutable.ListBuffer[StandardPersistencePair]
    vertFiltration.foreach { n1 =>
      val neighbors = link(n1).toList
      neighbors.size match {
        case 0 => {
          componentRep(n1) = n1
          extrema += n1
        }
        case 1 => {
          val n2 = neighbors.head
          val rep = componentRep(uf.find(n2))
          uf.union(n1, n2)
          componentRep(uf.find(n1)) = rep
        }
        case _ => {
          for ((n2, n3) <- neighbors.zip(neighbors.drop(1))) {
            val p2 = uf.find(n2)
            val p3 = uf.find(n3)
            if (p2 != p3) {
              // Two components merge.  Output a persistence pair
              val c = comp(componentRep(p2), componentRep(p3))
              if (c < 0) {
                pairs += new StandardPersistencePair(componentRep(p3), n1, componentRep(p2))
                uf.union(n1, n2)
                uf.union(n1, n3)
                componentRep(uf.find(n1)) = componentRep(p2)
              } else if (c > 0) {
                pairs += new StandardPersistencePair(componentRep(p2), n1, componentRep(p3))
                uf.union(n1, n2)
                uf.union(n1, n3)
                componentRep(uf.find(n1)) = componentRep(p3)
              } else {
                throw new RuntimeException
              }
            } else {
              val rep = componentRep(uf.find(n2))
              uf.union(n1, n2)
              uf.union(n1, n3)
              componentRep(uf.find(n1)) = rep
            }
          }
        }
        //case n => {throw new RuntimeException("Too many lower neighbors during merge: " + n)}
      }
    }

    println("numExtrema:" + extrema.size)

    pairs.toList
  }
}