/**
 * Ayla Geometry Core
 * (c) 2011-2014 William Harvey
 * http://www.aylasoftware.org
 */
package ayla.geometry

abstract class StandardPersistencePairLike(val extremum: Int, val saddle: Int, val killedBy: Int) extends Ordered[StandardPersistencePairLike] {
	val persistence: Float
	override def compare(o: StandardPersistencePairLike) = persistence.compare(o.persistence)
}
