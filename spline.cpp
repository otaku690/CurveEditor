#include "spline.h"

const int Division = 20;

Spline::Spline(GraphWidget *graphWidget)
:graph(graphWidget)
{
	deleteNodes(interpPoints);
	deleteNodes(controlPoints);
	deleteNodes(curvePoints);
    
	deleteEdges(pseudoEdges);
	deleteEdges(curveEdges);
	deleteEdges(controlEdges);
    
	startPoint = NULL;
	endPoint = NULL;
	totalPoints = 0;
	style = BezierBernstein;
	showControlPoints = false;
	showCurvePoints = false; // LOOK
	showPseudoPoints = true;
}

void Spline::SetVisible(QList<Node *> nodes, bool visible)
{
	foreach(Node * i, nodes)
	{
		i->setVisible(visible);
		i->setEnabled(visible);
	}
}
void Spline::SetVisible(QList<Edge *> edges, bool visible)
{
	foreach(Edge * i, edges)
	{
		i->setVisible(visible);
		i->setEnabled(visible);
	}
}
void Spline::SetVisible(bool visible)
{
	if (this->startPoint && this->endPoint)
	{
		this->startPoint->setVisible(visible);
		this->startPoint->setEnabled(visible);
		this->endPoint->setVisible(visible);
		this->endPoint->setEnabled(visible);
		if (!visible)
		{
			this->graph->scene()->removeItem(startPoint);
			this->graph->scene()->removeItem(endPoint);
		}
		else
		{
			this->graph->scene()->addItem(startPoint);
			this->graph->scene()->addItem(endPoint);
		}
	}
}

void Spline::Clear()
{
	deleteNodes(interpPoints);
	deleteNodes(controlPoints);
	deleteNodes(curvePoints);
    
	deleteEdges(pseudoEdges);
	deleteEdges(curveEdges);
	deleteEdges(controlEdges);
    
	if (startPoint)
	{
		this->graph->scene()->removeItem(startPoint);
		delete startPoint;
		startPoint = NULL;
	}
	if (endPoint)
	{
		this->graph->scene()->removeItem(endPoint);
		delete endPoint;
		endPoint = NULL;
	}
    
	totalPoints = 0;
}

void Spline::PaintPseudeLines()
{
	if (startPoint && endPoint)
	{
		Interpolate();
		deleteEdges(pseudoEdges);
		calcPseudoEdges();
		Paint();
	}
}

void Spline::deleteNodes(QList<Node *> &nodes)
{
	foreach(Node * i, nodes)
	{
		delete i;
	}
	nodes.clear();
}

void Spline::deleteEdges(QList<Edge *> &edges)
{
	foreach(Edge * i, edges)
	{
		delete i;
	}
	edges.clear();
}

void Spline::Paint()
{
	deleteEdges(curveEdges);
	deleteEdges(pseudoEdges);
	deleteEdges(controlEdges);
    
	// Curve Points & Edges
    // LOOK
	if(showCurvePoints)
	{
		foreach(Node * i, curvePoints)
		{
			this->graph->scene()->addItem(i);
		}
	}
	calcCurveEdges();
    
	// Control Points & Edges
	if (showControlPoints)
	{
		foreach(Node * i, controlPoints)
		{
			this->graph->scene()->addItem(i);
		}
		calcControlEdges();
	}
    
	if (showPseudoPoints)
	{
		calcPseudoEdges();
		this->SetVisible(true);
	}
	else
	{
		this->SetVisible(false);
	}
    
	this->graph->repaint(-1000, -1000, 2000, 2000);
}

void Spline::AddPoint(Node *newNode)
{
	interpPoints.append(newNode);
	graph->scene()->addItem(newNode);
	totalPoints++;
    
	CalcStartPoint();
	CalcEndPoint();
}

// LOOK (Adding intermediate control points)
void Spline::AddPoint(Node *newNode, int& index)
{
	// Add point as above but using the index value
}

// LOOK (Adding intermediate control points)
void Spline::DeletePoint(int& index)
{
	// Remove Node
    
	// Remove the start and end points if points less than 2
	if (totalPoints<=2)
	{
		if (startPoint)
		{
			ClearStartPoint();
		}
		if (endPoint)
		{
			ClearEndPoint();
		}
	}
	else
	{
		// Calculate the first and last points if those get deleted
		if (index == 0)
		{
			ClearStartPoint();
			CalcStartPoint();
		}
		if(index==totalPoints)
		{
			ClearEndPoint();
			CalcEndPoint();
		}
	}
    
	// Recalculate start and end points if those are the only points
	if(totalPoints==2)
	{
		ClearStartPoint();
		ClearEndPoint();
		CalcStartPoint();
		CalcEndPoint();
	}
}

void Spline::ClearStartPoint()
{
	this->graph->scene()->removeItem(startPoint);
	delete startPoint;
	startPoint = NULL;
}

void Spline::ClearEndPoint()
{
	this->graph->scene()->removeItem(endPoint);
	delete endPoint;
	endPoint = NULL;
}

void Spline::CalcStartPoint()
{
	if (totalPoints == 2 || (totalPoints>2 && !startPoint))
	{
		Node* p1 = interpPoints[0];
		Node* p2 = interpPoints[1];
        
		QPointF line12 = p1->pos() - p2->pos();
		double length12 = sqrt(line12.x()*line12.x() + line12.y()*line12.y());
        
		startPoint = new Node(this->graph, PseudoPoints);
		startPoint->setPos(p1->pos() + line12 / length12 * 50);
		graph->scene()->addItem(startPoint);
	}
}

void Spline::CalcEndPoint()
{
	if (totalPoints < 2)
		return;
	if (!endPoint)
	{
		endPoint = new Node(this->graph, PseudoPoints);
		graph->scene()->addItem(endPoint);
	}
    
	Node* pN = interpPoints[totalPoints - 1];
	Node* pNm1 = interpPoints[totalPoints - 2];
    
	QPointF lineNNm1 = pN->pos() - pNm1->pos();
	double lengthNNm1 = sqrt(lineNNm1.x()*lineNNm1.x() + lineNNm1.y()*lineNNm1.y());
    
	endPoint->setPos(pN->pos() + lineNNm1 / lengthNNm1 * 50);
}

void Spline::calcCurveEdges()
{
	deleteEdges(curveEdges);
	for (int i = 0; i < curvePoints.size() - 1; i++)
	{
		Edge * tmp = new Edge(curvePoints[i]->pos(), curvePoints[i+1]->pos());
		curveEdges.append(tmp);
		graph->scene()->addItem(tmp);
	}
}

void Spline::calcControlEdges()
{
	deleteEdges(controlEdges);
    
	switch(style)
	{
        case BezierBernstein:
        case BezierDeCasteljau:
        case BezierMatrix:
            if (this->controlPoints.size() == 2*(totalPoints-1))
            {
                for (int i = 0; i < this->totalPoints - 1; i ++)
                {
                    Node * n1 = interpPoints[i];
                    Node * n2 = interpPoints[i+1];
                    Node * c1 = controlPoints[2*i];
                    Node * c2 = controlPoints[2*i+1];
                    Edge * e1 = new Edge(n1->pos(), c1->pos(), 1);
                    Edge * e2 = new Edge(c1->pos(), c2->pos(), 1);
                    Edge * e3 = new Edge(c2->pos(), n2->pos(), 1);
                    controlEdges.append(e1);
                    controlEdges.append(e2);
                    controlEdges.append(e3);
                    graph->scene()->addItem(e1);
                    graph->scene()->addItem(e2);
                    graph->scene()->addItem(e3);
                }
            }
            break;
        case BSpline:
            for (int i = 0; i < this->controlPoints.size() - 1; i++)
            {
                Node * c1 = controlPoints[i];
                Node * c2 = controlPoints[i+1];
                
                Edge * e1 = new Edge(c1->pos(), c2->pos(), 1);
                controlEdges.append(e1);
                graph->scene()->addItem(e1);
            }
            break;
        case HermiteClamped:
        case HermiteNatural:
            if (this->controlPoints.size() == totalPoints)
            {
                for (int i = 0; i < this->totalPoints; i++)
                {
                    Node * c1 = interpPoints[i];
                    Node * c2 = controlPoints[i];
                    
                    c2->setPosition(c1->getPosition() + c2->getPosition());
                    
                    Edge * e1 = new Edge(c1->pos(), c2->pos(), 1);
                    controlEdges.append(e1);
                    graph->scene()->addItem(e1);
                }
            }
            break;
	}
}

void Spline::calcPseudoEdges()
{
	if (totalPoints >= 2)
	{
		deleteEdges(pseudoEdges);
		// start
		Edge * start = new Edge(startPoint->pos(), ((Node *)interpPoints.first())->pos(), 2);
		pseudoEdges.append(start);
		graph->scene()->addItem(start);
		// end
		Edge * end = new Edge(endPoint->pos(), ((Node *)interpPoints.last())->pos(), 2);
		pseudoEdges.append(end);
		graph->scene()->addItem(end);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Start your code here.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Spline::Interpolate()
{
	//Clear the old curve points
	deleteNodes(curvePoints);
	//Calculate control points.
	CalcControlPoints();
	//Depending on the selected style, interpolate the curve.
	switch(style){
        case BezierBernstein:		InterpBernstein();       break;
        case BezierDeCasteljau:		InterpCasteljau();      break;
        case BezierMatrix:          InterpMatrix();         break;
        case BSpline:               InterpBSpline();        break;
        case HermiteClamped:        InterpHermiteClamped(); break;
        case HermiteNatural:        InterpHermiteNatural(); break;
	}
}


//////////////////////////////////////////////////////////////////////////
// Calculate the control points
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// controlPoints- type: QList<Node *>
//				  discription: stores all the control points for curve interpolation and display.
//                             For Bezier curve, between very pair of consecutive interpolation points,
//                             there should be two control points. These four points determins the curve interpolation.
//                             For B-Spline, there should be totalPoints + 2 control points calculated from Ac = p.
//
// Hint: To implement Hermite and B-Spline, you need to write functions to create the A matrix as in the handouts.
//       Then you solve a linear system Ac = p, where p is the interpolation points vector and c are the control points.
//       We have provided you with a data structure to store and solve the linear system.
//       Below is an example code, read the understand it.
//
//	matrix<float> A(3,3);
//  matrix<float> c(3,1);
//  matrix<float> p(3,1);
//  A(0,0) = 1.0; A(0,1) = 0.0; A(0,2) = 0.0;
//  A(1,0) = 0.0; A(1,1) = 1.0; A(1,2) = 0.0;
//  A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 1.0;
//  p(0,0) = 1.0; p(1,0) = 2.0; p(2,0) = 3.0;
//  c = A.Solve(p);
//
//  The result in c is c(0,0) = 1.0; c(1,0) = 2.0; c(3,0) = 3.0, which satisfies Ac = p.
// Hint2: To get the position which is a vec3 in the class Node*, we need to use the method Node::getPosition().
// Example:
// vec3 startPointPos = startPoint->getPosition();
//
// Hint3: To new a Node*, add it into the List, and set the position of it, you need to do like this:
// Example:
// Node* ctrl = new Node(graph, ControlPoints);
// ctrl->setPosition(/*some vec3*/);
// controlPoints.append(ctrl);
void Spline::CalcControlPoints()
{
    float ratio = 0.33333333333;
    Node* pCP[2];
    vec3 pos;
    Node* p1,* p2;
	// Clear previous controlPoints
	deleteNodes(controlPoints);
	// Based on style use appropriate ways to compute control points
	switch(style)
	{
        case BezierBernstein:
        case BezierDeCasteljau:
        case BezierMatrix:

        if( interpPoints.size() < 2 )
            break;

        for( int i = 0; i < interpPoints.size()-1; ++i )
        {
            for( int c = 0; c < 2; ++c )
            {
                p1 = i+c-1 < 0 ? startPoint : interpPoints[i+c-1];
                p2 = i+c+1 >= interpPoints.size() ? endPoint : interpPoints[i+c+1];
                pCP[c] = new Node( graph, ControlPoints );
                //                    pos[0] = interpPoints[i+c]->x() + ratio * 0.5 *
                //                                ( p2->x() + p1->x() );

                //                    pos[1] = interpPoints[i+c]->y() + ratio * 0.5 *
                //                            ( p2->y() + p1->y() );
                pos = interpPoints[i+c]->getPosition() + ratio * 0.5 * ( p2->getPosition() - p1->getPosition() );

                pCP[c]->setPosition( pos );
                controlPoints.append(pCP[c]);
                ratio = -ratio;
            }
        }
        break;

        case BSpline:
            // make sure there are at least 2 totalPoints.
            // ADD YOUR CODE HERE
            break;
            // In the case for Hermite, you want to implement both "clamped" and "natural" versions.
        case HermiteClamped:
            // make sure there are at least 2 totalPoints.
            // ADD YOUR CODE HERE
            break;
        case HermiteNatural:
            // make sure there are at least 2 totalPoints.
            // ADD YOUR CODE HERE
            break;
	}
}

inline int factorial( int a )
{
    if( a < 2 ) return 1;
    int total =1;
    for(int i=a; i>1;--i) total *= i;
    return total;
}


//////////////////////////////////////////////////////////////////////////
// Cubic Berstein Bezier Spline
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points
//
// Hint: To get the position which is a vec3 in the class Node*, we need to use the method Node::getPosition().
// Example:
// vec3 b = interpPoints[i]->getPosition();
//
// Hint2: To new a Node*, add it into the List, and set the position of it, you need to do like this:
// Example:
// Node* newNode = new Node(graph, CurvePoints);
// newNode->setPosition(/*some vec3*/);
// curvePoints.append(newNode);

void Spline::InterpBernstein()
{
    if( controlPoints.size() < 2 )
        return;
    float delta = 1.0f/Division;
    float uu;
    Node* cP[4];
    Node* curvePoint;
    vec3 pos;
    for( int i = 0; i < interpPoints.size()-1; ++i )
    {
        cP[0] = interpPoints[i];
        cP[1] = controlPoints[2*i];
        cP[2] = controlPoints[2*i+1];
        cP[3] = interpPoints[i+1];

        uu = delta;
        for( int u = 1; u < Division; ++u )
        {
            curvePoint = new Node( graph, CurvePoints );

            pos = vec3(0);
            for( int deg= 0; deg <= 3; ++deg )
            {
                pos += pow( uu, deg ) * pow( 1-uu, 3-deg ) *
                        6.0f / ( factorial(deg) * factorial( 3-deg ) ) * cP[deg]->getPosition();

            }

            curvePoint->setPosition( pos );
            curvePoints.append( curvePoint );
            uu += delta;
        }

    }

}

//////////////////////////////////////////////////////////////////////////
// Cubic de Casteljau Bezier Spline
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points
//
// Hint: To get the position which is a vec3 in the class Node*, we need to use the method Node::getPosition().
// Example:
// vec3 b = interpPoints[i]->getPosition();
//
// Hint2: To new a Node*, add it into the List, and set the position of it, you need to do like this:
// Example:
// Node* newNode = new Node(graph, CurvePoints);
// newNode->setPosition(/*some vec3*/);
// curvePoints.append(newNode);
void Spline::InterpCasteljau()
{
	// ADD YOU CODE HERE
}

//////////////////////////////////////////////////////////////////////////
// Interpolation Matrix Bezier Spline
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points
//
// Hint: To get the position which is a vec3 in the class Node*, we need to use the method Node::getPosition().
// Example:
// vec3 b = interpPoints[i]->getPosition();
//
// Hint2: To new a Node*, add it into the List, and set the position of it, you need to do like this:
// Example:
// Node* newNode = new Node(graph, CurvePoints);
// newNode->setPosition(/*some vec3*/);
// curvePoints.append(newNode);
// Hint3: If you have any problem in vec3 vec4 mat4 operations, see the SvzAlgebra.h and svzalgebra.cpp to find some specific method out.
// For Example:
// vec3 operator * (const mat4& a, const vec3& v); has been defined in SvzAlgebra.h
// That means you can do operations like this:
// vec3 A;
// mat4 M;
// vec3 B = M*A;
// (It automatically transfer the A and B vectors to homogeneous vectors)
float m[16] = { -1.0f, 3.0f, -3.0f, 1.0f,
                 3.0f, -6.0f, 3.0f, 0,
                -3.0f, 3.0f, 0.0f, 0,
                 1.0f, 0, 0, 0 };


void Spline::InterpMatrix()
{
    if( controlPoints.size() < 2 )
        return;
    float delta = 1.0f/Division;
    float uu;
    vec3 cP[4];
    Node* curvePoint;
    vec3 pos;
    vec4 uvec;
    mat4 m4( m );

    for( int i = 0; i < interpPoints.size()-1; ++i )
    {
        cP[0] = interpPoints[i]->getPosition();
        cP[1] = controlPoints[2*i]->getPosition();
        cP[2] = controlPoints[2*i+1]->getPosition();
        cP[3] = interpPoints[i+1]->getPosition();

        uu = delta;
        for( int u = 1; u < Division; ++u )
        {
            curvePoint = new Node( graph, CurvePoints );
            uvec = vec4( uu*uu*uu, uu*uu, uu, 1 );
            pos[0] = uvec * m4 * vec4( cP[0][0], cP[1][0], cP[2][0], cP[3][0] );
            pos[1] = uvec * m4 * vec4( cP[0][1], cP[1][1], cP[2][1], cP[3][1] );
            pos[2] = uvec * m4 * vec4( cP[0][2], cP[1][2], cP[2][2], cP[3][2] );


            curvePoint->setPosition( pos );
            curvePoints.append( curvePoint );
            uu += delta;
        }

    }

}

//////////////////////////////////////////////////////////////////////////
// BSpline curve
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points
void Spline::InterpBSpline()
{
	// ADD YOUR CODE HERE
}

//////////////////////////////////////////////////////////////////////////
// Hermite Curve ('Clamped' condition)
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points
void Spline::InterpHermiteClamped()
{
	// ADD YOUR CODE HERE
}

//////////////////////////////////////////////////////////////////////////
// Hermite Curve ('Natural' condition)
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points 
// Hint: It is exactly the same like the method above.
void Spline::InterpHermiteNatural()
{
	// ADD YOUR CODE HERE.
}
