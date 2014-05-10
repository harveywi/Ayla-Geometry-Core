/**
 * Ayla Geometry Core
 * (c) 2011-2014 William Harvey
 * http://www.aylasoftware.org
 */
package ayla.geometry.landscape

import java.awt.geom.PathIterator
import java.awt.geom.Line2D
import java.awt.geom.Area
import java.awt.geom.Path2D
import java.awt.Polygon
import ayla.geometry.contourtree._
import ayla.geometry._
import scala.collection.mutable.ArrayBuffer
import java.awt.geom.Rectangle2D
import java.awt.geom.GeneralPath
import java.awt.geom.AffineTransform
import java.awt.image.BufferedImage
import java.awt.Graphics2D
import java.awt.Color
import java.awt.geom.AffineTransform
import org.jgrapht.graph._
import org.jgrapht.alg._
import scala.collection.JavaConversions._
import org.jgrapht.graph.{ DefaultDirectedGraph, DefaultEdge }
import scala.collection.mutable.Stack
import java.awt.geom.Point2D
import scala.collection.mutable.MultiMap
import ayla.util.Timing

object VoronoiFloorplan {
  val terrainHeight = .7f
  
  def apply(
    bspEnsemble: BSPTreeEnsemble,
    ct: ContourTree,
    epsilon: Double = 1e-4,
    delta: Double = 1e-5,
    pctChange: Double = 1d): ScalarFunction = {
    
    val outerBdry = Array((0, 0), (1, 0), (1, 1), (0, 1)).map { case (x, y) => new Point2D.Double(x, y).asInstanceOf[Point2D] }
    
    class VoronoiSplitRule extends SplitRule {
      def getNextSplit = this
    }
    
    val funcVals = ct.nodesContracted.map(n => ct.scalarFunction.getFuncVal(n.vertex))
    val minFuncVal = funcVals.min
    val maxFuncVal = funcVals.max
    val rangeFuncVal: Float = if (maxFuncVal == minFuncVal) 1 else maxFuncVal - minFuncVal
    
    def getVertexHeight(vertIdx: Int): Double = {
      val funcVal = ct.scalarFunction.getFuncVal(vertIdx)
      val height = terrainHeight * (funcVal - minFuncVal) / rangeFuncVal
      height
    }
    
    val floorplan = Floorplan[Array[Point2D]](bspEnsemble, scaleContour,
      splitContourCirclePartition, outerBdry, new VoronoiSplitRule, _ => "")
      
    buildTriangleArray(floorplan, getVertexHeight)
  }

  def buildTriangleArray(floorplan: Floorplan[Array[Point2D]], getVertexHeight: Int => Double): ScalarFunction = {

    val verts = new ArrayBuffer[Array[Float]]
    val stack = new Stack[FloorplanElem[Array[Point2D]]]

    val terrainVertToContourTreeNode = new ArrayBuffer[Int]

    class Subset(val range: Range, val e: ContourTreeEdge)
    val subsets = new scala.collection.mutable.ArrayBuffer[Subset]

    stack.pushAll(floorplan.rootElems)
    while (!stack.isEmpty) {
      val floorplanElem = stack.pop
      val contour = floorplanElem.contour
      if (contour.isEmpty) {
        stack.pushAll(floorplanElem.children)
      } else {
        val hOuter = getVertexHeight(floorplanElem.bspNode.outerNode.vertex)
        val hInner = getVertexHeight(floorplanElem.bspNode.innerNode.vertex)
        val numNewVerts = if (floorplanElem.children.isEmpty) {
          // Triangulate terminal region
          triangulateTerminalPolygon(contour, verts, hOuter, hInner, floorplanElem.bspNode.ctEdge.noncriticalNodes, getVertexHeight, terrainVertToContourTreeNode,
            floorplanElem.bspNode.outerNode.vertex, floorplanElem.bspNode.innerNode.vertex)
        } else {
          // Triangulate an annulus
          val sumChildAreas = floorplanElem.bspNode.children.map(_.sumArea).sum
          val scaleFactor = math.sqrt(sumChildAreas) / math.sqrt(floorplanElem.bspNode.sumArea)
          val innerContour = scaleContour(contour, scaleFactor)

          stack.pushAll(floorplanElem.children)
          triangulateAnnulus(contour, innerContour, verts, hOuter, hInner, floorplanElem.bspNode.ctEdge.noncriticalNodes, getVertexHeight, terrainVertToContourTreeNode,
            floorplanElem.bspNode.outerNode.vertex, floorplanElem.bspNode.innerNode.vertex)
        }

        val subset = new Subset(Range(verts.size - numNewVerts, verts.size), floorplanElem.bspNode.ctEdge)
        subsets += subset
      }
    }

    val vertsDistinct = verts.map { case Array(x, y, z) => (x, y, z) }.distinct

    val vertToIndex = vertsDistinct.zipWithIndex.toMap

    val faces = verts.grouped(3).map { tri =>
    	tri.map(v => vertToIndex((v(0), v(1), v(2)))).toArray
    }.toArray
    
    val sc = new SimplicialComplex(vertsDistinct.map{case (v1, v2, v3) => Array[Float](v1, v2, v3)}.toArray, faces)
    
    val funcVals = vertsDistinct.map(_._3).toArray
    
    new ScalarFunction(sc, funcVals)
  }

  def triangulateTerminalPolygon(Pinit: Array[Point2D], quadVertices: ArrayBuffer[Array[Float]], zOuter: Double, zInner: Double,
    noncriticalNodes: List[ContourTreeNode], getVertexHeight: Int => Double, terrainVertToContourTreeNode: ArrayBuffer[Int],
    outerNode: Int, innerNode: Int): Int = {
		val scaleFactors = Array(.95, .90, .80, .70, .55, .20, .05)
		val alphaVals = Array(.95, .85, .65, .4, .2, .1, .03)
//    val scaleFactors = Array(.5)
//    val alphaVals = Array(.5)

    val numLayers = scaleFactors.size

    var alphaOut = 1d
    var P = Pinit
    var totalVertsAdded = 0
    scaleFactors.zip(alphaVals).foreach {
      case (scaleFactor, alphaIn) => {
        val Q = scaleContour(P, scaleFactor)

        val zOut = alphaOut * zOuter + (1 - alphaOut) * zInner
        val zIn = alphaIn * zOuter + (1 - alphaIn) * zInner

        val zMin = math.min(zOut, zIn)
        val zMax = math.max(zOut, zIn)

        val ncStratum = noncriticalNodes.flatMap(n => {
          val z = getVertexHeight(n.vertex)
          if (z >= zMin && z < zMax)
            Some(n)
          else
            None
        })

        val nOuter = if (scaleFactor == .95) outerNode else -1

        totalVertsAdded += triangulateAnnulus(P, Q, quadVertices, zOut, zIn, ncStratum, getVertexHeight, terrainVertToContourTreeNode, nOuter, -1)
        alphaOut = alphaIn
        P = Q
      }
    }

    def makePoint3f(xy: Point2D, d3: Double) = {
      assert(xy.getX == xy.getX)
      assert(xy.getY == xy.getY)
      assert(d3 == d3)
      Array[Float](xy.getX.toFloat - .5f, xy.getY.toFloat - .5f, d3.toFloat)
    }

    val center = PolygonUtilities.centerOfMass(P)
    if (center.getX.isInfinity || center.getY.isInfinity) {
      val xSum = P.map(_.getX).sum
      val ySum = P.map(_.getY).sum
      center.setLocation(xSum / P.size.toDouble, ySum / P.size.toDouble)
    }
    val zOut = alphaOut * zOuter + (1 - alphaOut) * zInner
    (0 until P.size).foreach(i => {
      val j = (i + 1) % P.size
      val p1 = makePoint3f(P(i), zOut)
      val p2 = makePoint3f(P(j), zOut)
      val p3 = makePoint3f(center, zInner)
      quadVertices ++= List(p1, p2, p3)
      terrainVertToContourTreeNode ++= List(-1, -1, innerNode)
      totalVertsAdded += 3
    })

    // Stuff for computing the spheres
    val zMin = math.min(zOut, zInner)
    val zMax = math.max(zOut, zInner)

    val ncStratum = noncriticalNodes.flatMap(n => {
      val z = getVertexHeight(n.vertex)
      if (z >= zMin && z < zMax)
        Some(n)
      else
        None
    })

    return totalVertsAdded
  }

  def triangulateAnnulus(P: Array[Point2D], Q: Array[Point2D], quadVertices: ArrayBuffer[Array[Float]], zOuter: Double,
    zInner: Double, noncriticalNodes: List[ContourTreeNode], getVertexHeight: Int => Double, terrainVertToContourTreeNode: ArrayBuffer[Int],
    outerNode: Int, innerNode: Int): Int = {
    assert(P.size == Q.size)

    def makePoint3f(xy: Point2D, d3: Double) = {
      assert(xy.getX == xy.getX)
      assert(xy.getY == xy.getY)
      assert(d3 == d3)
      
      Array[Float](xy.getX.toFloat - .5f, xy.getY.toFloat - .5f, d3.toFloat)
    } 

    var totalVertsAdded = 0

    (0 until P.size).foreach(i => {
      val j = (i + 1) % P.size

      val ul = makePoint3f(P(i), zOuter)
      val ur = makePoint3f(P(j), zOuter)
      val ll = makePoint3f(Q(i), zInner)
      val lr = makePoint3f(Q(j), zInner)

      // Counter-clockwise
      quadVertices ++= List(ur, ul, ll)
      terrainVertToContourTreeNode ++= List(outerNode, outerNode, innerNode)

      quadVertices ++= List(ll, lr, ur)
      terrainVertToContourTreeNode ++= List(innerNode, innerNode, outerNode)

      totalVertsAdded += 6
    })

    return totalVertsAdded
  }

  def scaleContour(contour: Array[Point2D], scaleFactor: Double): Array[Point2D] = {
    if (contour.isEmpty)
      return Array.empty[Point2D]
    val oldCentroid = PolygonUtilities.centerOfMass(contour)
    val contour2 = new Array[Point2D](contour.size)
    AffineTransform.getScaleInstance(scaleFactor, scaleFactor).transform(contour, 0, contour2, 0, contour.size)
    val newCentroid = PolygonUtilities.centerOfMass(contour2)

    val tx = oldCentroid.getX - newCentroid.getX
    val ty = oldCentroid.getY - newCentroid.getY
    AffineTransform.getTranslateInstance(tx, ty).transform(contour2, 0, contour2, 0, contour2.size)

    return contour2
  }

  def splitContourCirclePartition(bdry: Array[Point2D], alphas: Seq[Double], rule: SplitRule): Seq[Array[Point2D]] = {

    if (bdry.isEmpty)
      return Seq.fill(alphas.size)(Array.empty[Point2D])

    // Find the farthest pair of points
    class DistTriple(val i: Int, val j: Int, val dSq: Double) extends Ordered[DistTriple] {
      def compare(that: DistTriple): Int = dSq.compareTo(that.dSq)
    }
    val farthestPair = (0 until bdry.size).flatMap(i => {
      (i + 1 until bdry.size).map(j => {
        new DistTriple(i, j, bdry(i).distanceSq(bdry(j)))
      })
    }).max

    // Build vector from one pole to the other, and compute projections of all points onto this vector.
    val polePt1 = bdry(farthestPair.i)
    val polePt2 = bdry(farthestPair.j)
    val vPole = Array[Double](polePt2.getX - polePt1.getX,
      polePt2.getY - polePt1.getY)
    val vPoleNorm = math.sqrt(vPole(0) * vPole(0) + vPole(1) * vPole(1))
    vPole.indices.foreach(vPole(_) /= vPoleNorm)
    val vTheta = math.atan2(vPole(1), vPole(0))

    val trans = new AffineTransform
    trans.rotate(-vTheta)
    trans.translate(-polePt1.getX, -polePt1.getY)
    // Create a rotated version of the points
    val bdryRotated = bdry.map(p =>
      trans.transform(p, null) //			new Point2D.Double((p.getX - polePt1.getX)*math.cos(-vTheta), (p.getY - polePt1.getY)*math.sin(-vTheta))
      ).sortWith((p1, p2) => p1.getX < p2.getX)

    // Get each point onto the positive y-axis, adding in its opposite point
    val ptsYAxis = bdryRotated.zipWithIndex.map {
      case (p, i) =>
        if (i == 0) {
          p
        } else if (i == bdryRotated.size - 1) {
          p
        } else {
          val s = math.signum(p.getY)
          val pLeft = (i - 1 to 0 by -1).view.map(bdryRotated(_)).find(p => math.signum(p.getY) != s) match {
            case Some(point) => point
            case None => bdryRotated(0)
          }
          val pRight = (i + 1 until bdryRotated.size).view.map(bdryRotated(_)).find(p => math.signum(p.getY) != s) match {
            case Some(point) => point
            case None => bdryRotated(bdryRotated.size - 1)
          }

          // Calculate the new y-value
          val alpha = if (pRight.getX != pLeft.getX) (p.getX - pLeft.getX) / (pRight.getX - pLeft.getX) else 0d
          val yAdd = pLeft.getY + alpha * (pRight.getY - pLeft.getY)
          val yNew = math.abs(p.getY) + math.abs(yAdd)
          new Point2D.Double(p.getX, yNew)
        }
    }

    val trueArea = PolygonUtilities.getArea(bdry)
    if (trueArea.isNaN) {
      //			System.err.println("Warning:  Zero-area polygon")
      //			return Seq.empty[Array[Point2D]]
      return Seq.fill(alphas.size)(Array.empty[Point2D])
    }

    var slicedBdry = bdry

    val subContours = alphas.dropRight(1).map(alpha => {
      val targetArea = trueArea * (alpha)
      //			println("Alpha is " + alpha)
      // now find the split point
      var areaSum = 0d
      var xSplit = -1d
      (0 until bdry.size - 1).foreach(i => {
        val p1 = ptsYAxis(i)
        val p2 = ptsYAxis(i + 1)
        val dx = p2.getX - p1.getX
        val areaRight = areaSum + dx * (p2.getY + p1.getY) / 2.0
        if (targetArea >= areaSum && targetArea <= areaRight) {
          // Calculate the split point here
          val s = targetArea - areaSum
          val a = p1.getY
          val b = p2.getY
          val h = p2.getX - p1.getX
          def quadForm(a: Double, b: Double, c: Double) = (-b + math.sqrt(b * b - 4 * a * c)) / (2.0 * a)

          val xDelta = if (a == b) 0 else quadForm((b - a) / (2 * h), a, -s)
          xSplit = p1.getX + xDelta
          if (xSplit.isNaN) {
            val t = 3
          }
          val areaDebug = areaSum + xDelta * (p2.getY + p1.getY) / 2.0
          //					println("Desired area:  " + targetArea + "\tDebug area:  " + areaDebug)
        }
        areaSum = areaRight
      })

      assert(areaSum >= 0)

      val split1Pre = new Point2D.Double(xSplit, -10)
      val split2Pre = new Point2D.Double(xSplit, 10)
      val split1Post = trans.inverseTransform(split1Pre, null)
      val split2Post = trans.inverseTransform(split2Pre, null)

      splitPolygonAtLine(split1Post, split2Post, slicedBdry) match {
        case Some((c1, c2)) => {
          val area1 = PolygonUtilities.getArea(c1)
          val area2 = PolygonUtilities.getArea(c2)

          if (area1.isNaN || area2.isNaN) {
            val t = 3
          }
          //					println("\tPoly areas:  " + area1 + "\t" + area2)
          val frac1 = area1 / trueArea
          val frac2 = area2 / trueArea
          val alphaBdry = 1 - alpha
          //					println(List(frac1, frac2, alphaBdry).mkString(","))
          if (math.abs(alphaBdry - frac1) < math.abs(alphaBdry - frac2)) {
            slicedBdry = c1
            //						println("A")
            c2
          } else {
            //						println("B")
            slicedBdry = c2
            c1
          }
        }
        case None => {
          System.err.println("Warning:  Impossible split")
          Array.empty[Point2D]
        }
      }
    })
    if (slicedBdry.flatMap(p => List(p.getX, p.getY)).exists(_.isNaN)) {
      assert(false)
    }

    subContours.foreach(c => if (c.flatMap(p => List(p.getX, p.getY)).exists(_.isNaN)) {
      assert(false)
    })

    if (alphas.size > 2) {
      val z = 3
    }

    return subContours.toSeq ++ List(slicedBdry)
  }

  def splitContour(bdry: Array[Point2D], alpha: Double, rule: SplitRule, epsilon: Double, delta: Double, pctChange: Double): (Array[Point2D], Array[Point2D]) = {
    val (pct1, pct2, swap) = if (alpha < .5) {
      (alpha, 1 - alpha, false)
    } else {
      (1 - alpha, alpha, true)
    }

    val gp = getGeneralPath(bdry)
    val brect = gp.getBounds2D

    val bdryArea = PolygonUtilities.getArea(bdry)

    var p1 = getRandomPointInPolygon(gp, brect)
    var p2 = getRandomPointInPolygon(gp, brect)

    var w1 = .5d
    var w2 = .5d

    Stream.from(0).foreach(_ => {
      val sumWeight = w1 + w2
      val d = p1.distance(p2)
      if (sumWeight > d) {
        w1 = d * w1 / sumWeight
        w2 = d * w2 / sumWeight
      }

      var stable = true
      getVoronoiPolygons(bdry, p1, p2, w1, w2) match {
        case Some((v1, v2)) => {
          val areaCell1 = PolygonUtilities.getArea(v1)
          val areaCell2 = PolygonUtilities.getArea(v2)
          val pctAreaCell1 = areaCell1 / bdryArea
          val pctAreaCell2 = areaCell2 / bdryArea

          if (math.abs(pctAreaCell1 - pct1) > epsilon)
            stable = false

          // Adjust weights
          w1 += w1 * (pct1 - pctAreaCell1) / pct1 * pctChange
          w1 = math.max(w1, delta)

          w2 += w2 * (pct2 - pctAreaCell2) / pct2 * pctChange
          w2 = math.max(w2, delta)

          // Move generators
          val center1 = PolygonUtilities.centerOfMass(v1)
          val center2 = PolygonUtilities.centerOfMass(v2)

          var success = true
          success &= moveGenerator(gp, center1, p1, bdry)
          success &= moveGenerator(gp, center2, p2, bdry)

          if (!success) {
            println("moveGenerator failed.")
            p1 = getRandomPointInPolygon(gp, brect)
            p2 = getRandomPointInPolygon(gp, brect)
            val d = p1.distance(p2)
            val randPct1 = math.random
            val randPct2 = 1 - randPct1
            w1 = randPct1 * d
            w2 = randPct2 * d
            println("Resetting weights and generator locations...")
            stable = false
            // continue
          }
          if (stable) {
            println("Converged:  " + System.currentTimeMillis)
            if (swap) {
              return (v2, v1)
            } else {
              return (v1, v2)
            }
          }
        }
        case None => {
          p1 = getRandomPointInPolygon(gp, brect)
          p2 = getRandomPointInPolygon(gp, brect)
          val d = p1.distance(p2)
          val randPct1 = math.random
          val randPct2 = 1 - randPct1
          w1 = randPct1 * d
          w2 = randPct2 * d
          println("Resetting weights and generator locations...")
          stable = false
          // continue
        }
      }
    })

    null
  }

  def moveGenerator(gp: GeneralPath, center: Point2D, p: Point2D, bdry: Array[Point2D]): Boolean = {
    if (gp.contains(center)) {
      p.setLocation(center.getX, center.getY)
      return true
    } else {
      var nearestU = Double.MaxValue
      var nearestPt: Point2D = null
      val uOut = new Array[Double](2)
      (0 until bdry.length).foreach(i => {
        val j = (i + 1) % bdry.length
        val b1 = bdry(i)
        val b2 = bdry(j)
        getIntersectU(p, center, b1, b2, uOut) match {
          case Some(isect) => {
            if (uOut(0) < nearestU) {
              nearestU = uOut(0)
              nearestPt = isect
            }
          }
          case None => { /*Do nothing*/ }
        }
      })
      if (nearestPt == null)
        return false
      p.setLocation(nearestPt.getX, nearestPt.getY)
      return true
    }
  }
  def splitPolygonAtLine(lineEnd1: Point2D, lineEnd2: Point2D, bdry: Array[Point2D]): Option[(Array[Point2D], Array[Point2D])] = {
    var isect1: Point2D = null
    var idxIsect1: Int = -1
    var isect2: Point2D = null
    var idxIsect2: Int = -1
    val uOut = new Array[Double](2)
    (0 until bdry.length).foreach(i => {
      val j = (i + 1) % bdry.length
      getIntersectU(lineEnd1, lineEnd2, bdry(i), bdry(j), uOut) match {
        case Some(isect) => {
          if (isect1 == null) {
            isect1 = isect
            idxIsect1 = i
          } else {
            isect2 = isect
            idxIsect2 = i
          }
        }
        case None => { /*Do nothing*/ }
      }
    })

    if (isect1 == null || isect2 == null)
      return None

    val halfBdry1 = new ArrayBuffer[Point2D]
    val halfBdry2 = new ArrayBuffer[Point2D]
    var activeHalfBdry = halfBdry1

    (0 until bdry.size).foreach(i => {
      activeHalfBdry += bdry(i).clone.asInstanceOf[Point2D]
      if (i == idxIsect1 && activeHalfBdry == halfBdry1) {
        halfBdry1 += isect1
        halfBdry1 += isect2
        activeHalfBdry = halfBdry2
      } else if (i == idxIsect2 && activeHalfBdry == halfBdry2) {
        halfBdry2 += isect2
        halfBdry2 += isect1
        activeHalfBdry = halfBdry1
      } else if (i == idxIsect1 && activeHalfBdry == halfBdry2) {
        halfBdry2 += isect1
        halfBdry2 += isect2
        activeHalfBdry = halfBdry1
      } else if (i == idxIsect2 && activeHalfBdry == halfBdry1) {
        halfBdry1 += isect2
        halfBdry1 += isect1
        activeHalfBdry = halfBdry2
      }
    })

    if (halfBdry1.size < 3 || halfBdry2.size < 3) {
      println("Empty polygons.  Resetting weights...")
      return None
    }

    val subPoly1 = halfBdry1.toArray
    val subPoly2 = halfBdry2.toArray

    return Some(subPoly1, subPoly2)
  }

  def splitPolygonAtLine(q: Point2D, v: Array[Double], bdry: Array[Point2D]): Option[(Array[Point2D], Array[Point2D])] = {
    val lineEnd1 = new Point2D.Double(q.getX + v(0) * 9999, q.getY + v(1) * 9999)
    val lineEnd2 = new Point2D.Double(q.getX + v(0) * -9999, q.getY + v(1) * -9999)

    var isect1: Point2D = null
    var idxIsect1: Int = -1
    var isect2: Point2D = null
    var idxIsect2: Int = -1
    val uOut = new Array[Double](2)
    (0 until bdry.length).foreach(i => {
      val j = (i + 1) % bdry.length
      getIntersectU(q, lineEnd1, bdry(i), bdry(j), uOut) match {
        case Some(isect) => {
          isect1 = isect
          idxIsect1 = i
        }
        case None => { /*Do nothing*/ }
      }

      getIntersectU(q, lineEnd2, bdry(i), bdry(j), uOut) match {
        case Some(isect) => {
          isect2 = isect
          idxIsect2 = i
        }
        case None => { /*Do nothing*/ }
      }
    })

    if (isect1 == null || isect2 == null)
      return None

    val halfBdry1 = new ArrayBuffer[Point2D]
    val halfBdry2 = new ArrayBuffer[Point2D]
    var activeHalfBdry = halfBdry1

    (0 until bdry.size).foreach(i => {
      activeHalfBdry += bdry(i).clone.asInstanceOf[Point2D]
      if (i == idxIsect1 && activeHalfBdry == halfBdry1) {
        halfBdry1 += isect1
        halfBdry1 += isect2
        activeHalfBdry = halfBdry2
      } else if (i == idxIsect2 && activeHalfBdry == halfBdry2) {
        halfBdry2 += isect2
        halfBdry2 += isect1
        activeHalfBdry = halfBdry1
      } else if (i == idxIsect1 && activeHalfBdry == halfBdry2) {
        halfBdry2 += isect1
        halfBdry2 += isect2
        activeHalfBdry = halfBdry1
      } else if (i == idxIsect2 && activeHalfBdry == halfBdry1) {
        halfBdry1 += isect2
        halfBdry1 += isect1
        activeHalfBdry = halfBdry2
      }
    })

    if (halfBdry1.size < 3 || halfBdry2.size < 3) {
      println("Empty polygons.  Resetting weights...")
      return None
    }

    val subPoly1 = halfBdry1.toArray
    val subPoly2 = halfBdry2.toArray

    return Some(subPoly1, subPoly2)
  }

  def getVoronoiPolygons(bdry: Array[Point2D], p1: Point2D, p2: Point2D, w1: Double, w2: Double): Option[(Array[Point2D], Array[Point2D])] = {
    // Compute bisector
    val d = p1.distance(p2)
    val t = (d + w1 - w2) / 2.0

    val v = Array(p2.getX - p1.getX, p2.getY - p1.getY)
    val lengthV = math.sqrt(v(0) * v(0) + v(1) * v(1))
    v(0) /= lengthV
    v(1) /= lengthV

    val q = new Point2D.Double(p1.getX + v(0) * t, p1.getY + v(1) * t)
    val vPerp = Array(-v(1), v(0))

    splitPolygonAtLine(q, vPerp, bdry) match {
      case Some((c1, c2)) => {
        val gp1 = getGeneralPath(c1)
        if (gp1.contains(p1))
          return Some(c1, c2)
        else
          return Some(c2, c1)
      }
      case None => return None
    }
  }

  def getIntersectU(p1: Point2D, p2: Point2D, p3: Point2D, p4: Point2D, uOut: Array[Double]): Option[Point2D] = {
    val denom = (p4.getY - p3.getY) * (p2.getX - p1.getX) - (p4.getX - p3.getX) * (p2.getY - p1.getY)

    if (denom == 0) {
      uOut(0) = 0
      uOut(1) = 0
      return None
    }

    val numeratorA = (p4.getX - p3.getX) * (p1.getY - p3.getY) - (p4.getY - p3.getY) * (p1.getX - p3.getX)
    val numeratorB = (p2.getX - p1.getX) * (p1.getY - p3.getY) - (p2.getY - p1.getY) * (p1.getX - p3.getX)

    val uA = numeratorA / denom
    val uB = numeratorB / denom

    uOut(0) = uA
    uOut(1) = uB
    if (uA < 0 || uA > 1 || uB < 0 || uB > 1)
      return None

    val x = p3.getX + uB * (p4.getX - p3.getX)
    val y = p3.getY + uB * (p4.getY - p3.getY)

    return Some(new Point2D.Double(x, y))
  }

  def getGeneralPath(poly: Array[Point2D]): GeneralPath = {
    val gp = new GeneralPath
    gp.moveTo(poly(0).getX, poly(0).getY)
    (1 until poly.size).foreach(i => {
      gp.lineTo(poly(i).getX, poly(i).getY)
    })
    gp.closePath
    return gp
  }

  def getRandomPointInPolygon(gp: GeneralPath, brect: Rectangle2D): Point2D = {
    val p = new Point2D.Double
    Stream.from(0).foreach(_ => {
      p.x = brect.getX + math.random * brect.getWidth
      p.y = brect.getY + math.random * brect.getHeight
      if (gp.contains(p))
        return p
    })
    null
  }
}
