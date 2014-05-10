/**
 * Ayla Geometry Core
 * (c) 2011-2014 William Harvey
 * http://www.aylasoftware.org
 */
package ayla.geometry.contourtree

class ContourTreeEdge(
  val n1: ContourTreeNode,
  val n2: ContourTreeNode,
  val noncriticalNodes: List[ContourTreeNode]) extends Serializable {

  var area: Double = 0
  override def toString = "(" + n1.toString + ", " + n2.toString + ")"
}