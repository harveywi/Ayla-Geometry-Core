/**
 * Ayla Geometry Core
 * (c) 2011-2014 William Harvey
 * http://www.aylasoftware.org
 */
package ayla.geometry.contourtree

import scala.collection.mutable.ArrayBuffer

class ContourTreeNode(val vertex: Int) {
  val parents = new ArrayBuffer[ContourTreeNode]
  val children = new ArrayBuffer[ContourTreeNode]
  def isMax = (parents.size == 0 && children.size == 1)
  def isMin = (children.size == 0 && parents.size == 1)
  def isBranchUp = (parents.size > 1 && children.size == 1)
  def isBranchDown = (children.size > 1 && parents.size == 1)
  def isOrdinary = (parents.size == 1 && children.size == 1)
  def isMultibranch = (parents.size > 1 && children.size > 1)
  def isSaddle = (isBranchUp || isBranchDown)
  def isCritical = !isOrdinary
  override def toString = vertex.toString
}