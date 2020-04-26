import kotlin.math.acos
import kotlin.math.max
import kotlin.math.pow
import kotlin.math.sqrt

class Hole {
    val vertices = mutableListOf<Vertex>()
    val edges = mutableListOf<Edge>()
}

data class Weight(val maxDihedralAngle: Float, val area: Float) : Comparable<Weight> {
    constructor(mesh: Mesh, v1:Vertex, v2: Vertex, v3:Vertex) : this(maxDihedralAngle(mesh, v1, v2, v3), area(v1, v2, v3))
    override fun compareTo(other: Weight): Int {
        return compareValuesBy(this, other, {it.maxDihedralAngle}, {it.area})
    }
    operator fun plus(increment: Weight): Weight {
        return Weight(maxDihedralAngle + increment.maxDihedralAngle, area + increment.area)
    }
    override fun toString(): String = "%.2f".format(area*1000000)
}

// Planar equation: Ax + By + Cz = D
class Plane(v1:Vertex, v2: Vertex, v3:Vertex) {
    val A: Float
    val B: Float
    val C: Float
    init {
        val vector12 = listOf(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z)
        val vector13 = listOf(v3.x - v1.x, v3.y - v1.y, v3.z - v1.z)

        // cross product of two vectors => A*i + B*j + C*k
        A = vector12[1] * vector13[2] - vector12[2] * vector13[1]
        B = vector12[1] * vector13[2] - vector12[2] * vector13[1]
        C = vector12[0] * vector13[1] - vector12[1] * vector13[0]
    }
}

fun area(v1: Vertex, v2: Vertex, v3: Vertex) : Float {
    val a = sqrt((v1.x - v2.x).pow(2) + (v1.y - v2.y).pow(2) + (v1.z - v2.z).pow(2))
    val b = sqrt((v1.x - v3.x).pow(2) + (v1.y - v3.y).pow(2) + (v1.z - v3.z).pow(2))
    val c = sqrt((v2.x - v3.x).pow(2) + (v2.y - v3.y).pow(2) + (v2.z - v3.z).pow(2))
    val p = (a+b+c)/2
    return sqrt(p*(p-a)*(p-b)*(p-c))
}

// Dihedral angle between two planes.
// p = [A, B, C] where Ax + By + Cz = D
fun dihedralAngle(p1: Plane, p2: Plane): Float {
    return acos((p1.A * p2.A + p1.B * p2.B + p1.C * p2.C) /
            (sqrt(p1.A * p1.A + p1.B * p1.B + p1.C * p1.C) * sqrt(p2.A * p2.A + p2.B * p2.B + p2.C * p2.C)))
}

fun maxDihedralAngle(mesh: Mesh, v1: Vertex, v2: Vertex, v3: Vertex) : Float {
    // maximum dihedral angle between (v1, v2, v3) and existing adjacent triangles
    // collect all adjacent triangles
    val adjacentTriangles = mutableSetOf<Int>()
    listOf(v1, v2, v3).forEach {
        adjacentTriangles.addAll(it.adjacentTriangles)
    }
    val plane = Plane(v1, v2, v3)
    var maxDihedral = Float.MIN_VALUE
    for (tri in adjacentTriangles) {
        val t = mesh.triangles[tri]
        val plane2 = Plane(mesh.vertices[t.v1], mesh.vertices[t.v2], mesh.vertices[t.v3])
        maxDihedral = max(maxDihedral, dihedralAngle(plane, plane2))
    }
    return maxDihedral
}

// An edge is a boundary edge if its vertices share exactly one triangle
fun isBoundary(v1: Vertex, v2: Vertex): Boolean {
    return v1.adjacentTriangles.intersect(v2.adjacentTriangles).size == 1
}

fun isChecked(holes: List<Hole>, edgeId: Int) : Boolean {
    return holes.any { hole -> hole.edges.any { it.id == edgeId } }
}

fun identifyHoles(mesh: Mesh) : List<Hole> {
    val res = mutableListOf<Hole>()
    for(edge in mesh.edges) {
        if(!isChecked(res, edge.id) && isBoundary(mesh.vertices[edge.v1], mesh.vertices[edge.v2])) {
            val hole = Hole()
            var currentVertex = mesh.vertices[edge.v1]
            var currentEdge = edge
            do {
                hole.vertices.add(currentVertex)
                hole.edges.add(currentEdge)
                for (it in currentVertex.adjacentEdges) {
                    if(it != currentEdge.id) {
                        val v1 = mesh.vertices[ mesh.edges[it].v1 ]
                        val v2 = mesh.vertices[ mesh.edges[it].v2 ]
                        if(isBoundary(v1, v2)) {
                            currentEdge = mesh.edges[it]
                            currentVertex = if (v1 == currentVertex) v2 else v1
                            break
                        }
                    }
                }
            } while(edge.id != currentEdge.id)
            res.add(hole)
        }
    }
    return res
}

fun minWeightTable(mesh: Mesh, hole: Hole): Pair<MutableList<MutableList<Weight>>, MutableList<MutableList<Int>>> {
    //for each (Vi,Vj) pair in (V0..Vn) find the Vx that gives minimum weight(Vi, Vj, Vx)
    val n = hole.vertices.size
    val minWeightTable = MutableList(n) { MutableList(n) { Weight(0F, 0F) } }
    val minWeightIndexTable = MutableList(n) { MutableList(n) { -1 } }
    for(i in 0..n-3) {
        minWeightTable[i][i+2] = Weight(mesh, hole.vertices[i], hole.vertices[i+1], hole.vertices[i+2])
    }
    for(j in 2 until n) {
        for(i in 0 until n-j) {
            var k = i + j
            var weightIK = Weight(Float.MAX_VALUE, Float.MAX_VALUE)
            var lambdaIK = -1
            for (m in i + 1 until k) {
                val weightIMK = minWeightTable[i][m] +  minWeightTable[m][k] + Weight(mesh, hole.vertices[i], hole.vertices[m], hole.vertices[k])
                if (weightIMK < weightIK) {
                    weightIK = weightIMK
                    lambdaIK = m
                }
            }
            minWeightTable[i][j] = weightIK
            minWeightIndexTable[i][j] = lambdaIK
        }
    }
    return Pair(minWeightTable, minWeightIndexTable)
}

fun fillHole(mesh : Mesh, hole: Hole, minWeightTable: List<List<Int>>, begin: Int, end: Int) {
    if(end - begin > 1) {
        var e = end
        //Vertex V for which the weight of triangle (begin, V, end) is minimum
        var minWeightVertex = minWeightTable[begin][e]
        while(minWeightVertex !in begin until end){
            e--
            minWeightVertex = minWeightTable[begin][e]
        }
        mesh.addTriangle(hole.vertices[begin].id, hole.vertices[minWeightVertex].id, hole.vertices[end].id)
        fillHole(mesh, hole, minWeightTable, begin, minWeightVertex)
        fillHole(mesh, hole, minWeightTable, minWeightVertex, end)
    }
}

fun main(args : Array<String>) {

    if(args.size != 2) {
        println("Usage: HoleFiller input_file output_file")
        return
    }

    val mesh = Mesh()
    mesh.loadOff(args[0])

    val holes = identifyHoles(mesh)
    println("found ${holes.size} holes")

    // triangulate holes
    for (hole in holes) {
        println("triangulating hole of size ${hole.vertices.size}...")
        val minWeightTable = minWeightTable(mesh, hole)
        fillHole(mesh, hole, minWeightTable.second, 0, hole.vertices.size - 1)
    }

    if(holes.isNotEmpty()) {
        mesh.exportOff(args[1])
    }
}