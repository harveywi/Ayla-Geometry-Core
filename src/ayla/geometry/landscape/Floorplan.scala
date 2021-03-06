/**
 * Ayla Geometry Core
 * (c) 2011-2014 William Harvey
 * http://www.aylasoftware.org
 */
package ayla.geometry.landscape

import scala.annotation.tailrec

abstract class SplitRule {
	def getNextSplit: SplitRule
}

@SerialVersionUID(1L)
class FloorplanElem[T](val bspNode: BSPNode, val contour: T, val splitRule: SplitRule) extends Serializable {
	val children = new scala.collection.mutable.ArrayBuffer[FloorplanElem[T]]
}

@SerialVersionUID(1L)
class Floorplan[T](val rootElems: List[FloorplanElem[T]]) extends Serializable

object Floorplan {
	def apply[T](
			bspEnsemble: BSPTreeEnsemble,
			scaleContour: (T, Double) => T,
			splitContour: (T, Seq[Double], SplitRule) => Seq[T],
			outerBdry: T,
			defaultSplitRule: SplitRule,
			str: T => String): Floorplan[T] = {
		
		// Get the list of initial floorplan contours corresponding to the BSP roots
		val totalRootArea = bspEnsemble.roots.map(_.sumArea).sum
		val initialBoundaries = if (bspEnsemble.roots.size == 1)
			Seq(outerBdry)
		else {
			val alphas = bspEnsemble.roots.map(_.sumArea).scanLeft(0d)(_ + _).drop(1).map(_ / totalRootArea)
			splitContour(outerBdry, alphas, defaultSplitRule)
		}
		
		val floorplanElems = bspEnsemble.roots.zip(initialBoundaries).map{case (root, bdry) => new FloorplanElem[T](root, bdry, defaultSplitRule.getNextSplit)}
		
		// Perform the recursive partioning of each of the root contours
		val stack = new scala.collection.mutable.Stack[FloorplanElem[T]]
		stack.pushAll(floorplanElems)
		while (!stack.isEmpty) {
			val elem = stack.pop
			val sumChildAreas = elem.bspNode.children.map(_.sumArea).sum
			val scaleFactor = math.sqrt(sumChildAreas) / math.sqrt(elem.bspNode.sumArea)
			val innerContour = scaleContour(elem.contour, scaleFactor)
			
			val alphas = elem.bspNode.children.map(_.sumArea).scanLeft(0d)(_ + _).drop(1).map(_ / sumChildAreas)
			val subContours = splitContour(innerContour, alphas, elem.splitRule)
			assert(subContours.size == elem.bspNode.children.size)
			elem.children ++= elem.bspNode.children.zip(subContours).map{case (node, contour) => new FloorplanElem[T](node, contour, elem.splitRule.getNextSplit)}
			stack.pushAll(elem.children)
		}
		
		// Make sure that the areas are correct
		floorplanElems.foreach(e => {
			val bspNodeArea = e.bspNode.sumArea
			println("bspNode area: " + bspNodeArea)
		})
		
		return new Floorplan[T](floorplanElems)
	}
}
