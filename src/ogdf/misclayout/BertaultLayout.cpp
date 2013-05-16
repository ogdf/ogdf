
/** \file
 * \brief Declaration of class BertaultLayout.
 * Computes a force directed layout (Bertault Layout) for preserving the planar embedding in the graph.
 * The algorithm is based on the paper
 * "A force-directed algorithm that preserves
 * edge-crossing properties" by Francois Bertault
 *
 * \author Smit Sanghavi;
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation;
 * see the file LICENSE.txt included in the packaging of this file
 * for details.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/


#include <ogdf/misclayout/BertaultLayout.h>
#include <ogdf/basic/Math.h>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <ogdf/fileformats/GraphIO.h>

namespace ogdf {

BertaultLayout::BertaultLayout()
{
	req_length=0;
	iter_no=0;
		impred=false;
}

BertaultLayout::BertaultLayout(double length, int number)
{
	req_length=length;
	iter_no=number;
		impred=false;
}

BertaultLayout::BertaultLayout(int number)
{
	req_length=0;
	iter_no=number;
		impred=false;
}

BertaultLayout::~BertaultLayout()
{
}

void BertaultLayout::call(GraphAttributes &AG)
{
	const Graph &G = AG.constGraph();
	if(G.numberOfNodes() == 0)
		return;
	if( (AG.attributes() & GraphAttributes::nodeGraphics) == 0 )
		return;
	if( (AG.attributes() & GraphAttributes::edgeGraphics) != 0 )
		AG.clearAllBends();
	if(iter_no==0)
		iter_no=G.numberOfNodes()*10;
	if(req_length==0)
	{
		edge e;
		forall_edges(e,G)
		{
			node a=e->source();
			node b=e->target();
			req_length+=sqrt((AG.x(a)-AG.x(b))*(AG.x(a)-AG.x(b))+(AG.y(a)-AG.y(b))*(AG.y(a)-AG.y(b)));
		}
		req_length=req_length/(G.numberOfEdges());
	}
	limit=4*req_length;					// can be changed... this value is taken in the research paper
	F_x.init(G);
	F_y.init(G);
	sect.init(G);

	//impred=true;

	if(impred)
		preprocess(AG);

	for(int k=0;k<iter_no;k++)
	{
		node v;
		forall_nodes(v,G)
		{
				//initialise everything
			F_x[v]=0;
			F_y[v]=0;
			sect[v].initialize();
	   /*     int i;
			for(i=1;i<9;i++)
			{
				cout << sect[v].R[i];
				cout << " ";
			}
			cout << "\n";*/

		}

		forall_nodes(v,G)
		{
			node j;
			//calculate total node-node repulsive force
			forall_nodes(j,G)
			{
				if(j!=v)
				f_Node_Repulsive(&v,&j,AG);
			}

			//calculate total node-node attractive force
			adjEntry adj;
			forall_adj(adj,v)
			{
				node ad=adj->twinNode();
				f_Node_Attractive(&v,&ad,AG);
			}
			//calculate total node-edge repulsive force
			edge e;
			forall_edges(e,G)
			{
				if(e->target()!=v&&e->source()!=v)
				{

					compute_I(&v,&e,AG);			//computes the projection

					if(i_On_Edge(&e,AG))		//computes if projection is on the edge
					{
						if((!impred)||surr(v->index(),e->index())==1)
						{
							 f_Edge(&v,&e,AG);
						}
						r_Calc_On_Edge(&v,&e,AG);					// updates values of section radii
					}
					else
					{
						r_Calc_Outside_Edge(&v,&e,AG);				// updates values of section radii
					}

				}
			}
		}

		forall_nodes(v,G)
		{
			/*   cout << "F_x : ";
			cout << F_x[v];

			cout << "   F_y : ";
			cout << F_y[v];
			cout << "\n";
			*/
			//moves the nodes according to forces
			move(&v,AG);
			/*
			cout << "move_x : ";
			cout << F_x[v];

			cout << "   move_y : ";
			cout << F_y[v];
			cout << "\n";

			int i;
			for(i=1;i<9;i++)
			{
				cout << sect[v].R[i];
				cout << " ";
			}
			cout << "\n\n";*/
		}
	}

}


void BertaultLayout::f_Node_Repulsive(node *v,node *j, GraphAttributes &AG)
{
	double dist=sqrt((AG.x(*v)-AG.x(*j))*(AG.x(*v)-AG.x(*j))+(AG.y(*v)-AG.y(*j))*(AG.y(*v)-AG.y(*j)));
	(F_x)[*v]+=((req_length)/dist)*((req_length)/dist)*(AG.x(*v)-AG.x(*j));
	(F_y)[*v]+=((req_length)/dist)*((req_length)/dist)*(AG.y(*v)-AG.y(*j));
}

void BertaultLayout::f_Node_Attractive(node *v,node *j, GraphAttributes &AG)
{
	double dist=sqrt((AG.x(*v)-AG.x(*j))*(AG.x(*v)-AG.x(*j))+(AG.y(*v)-AG.y(*j))*(AG.y(*v)-AG.y(*j)));
	(F_x)[*v]+=(-(dist/req_length)*(AG.x(*v)-AG.x(*j)));
	(F_y)[*v]+=(-(dist/req_length)*(AG.y(*v)-AG.y(*j)));
}

void BertaultLayout::compute_I(node *v,edge *e, GraphAttributes &AG)
{
	node a=(*e)->source();
	node b=(*e)->target();
	double m=(AG.y(a)-AG.y(b))/(AG.x(a)-AG.x(b));			//slope of edge
	double n=-1/m;										//slope of a perpendicular
	double c=AG.y(a)-m*(AG.x(a));							//y=mx+c for edge
	double d=AG.y(*v)-n*(AG.x(*v));							//y=nx+d for the perpendicular
	i.x=(d-c)/(m-n);										//solve for x
	i.y=m*i.x+c;											//solve for y
}

bool BertaultLayout::i_On_Edge(edge *e, GraphAttributes &AG)
{
	node a=(*e)->source();
	node b=(*e)->target();
	return ((i.x<=AG.x(a)&&i.x>=AG.x(b))||(i.x>=AG.x(a)&&i.x<=AG.x(b)))&&((i.y<=AG.y(a)&&i.y>=AG.y(b))||(i.y>=AG.y(a)&&i.y<=AG.y(b)));				// x and y coordinates of i must be in between that of a and b
}

void BertaultLayout::f_Edge(node *v,edge *e, GraphAttributes &AG)
{
	double dist=sqrt((AG.x(*v)-i.x)*(AG.x(*v)-i.x)+(AG.y(*v)-i.y)*(AG.y(*v)-i.y));
	if(dist<=limit&&dist>0)
	{
		double fx=(limit-dist)*(limit-dist)*(AG.x(*v)-i.x)/dist;
		double fy=(limit-dist)*(limit-dist)*(AG.y(*v)-i.y)/dist;
		(F_x)[*v]+=fx;
		(F_y)[*v]+=fy;
		node a=(*e)->source();
		node b=(*e)->target();
		(F_x)[a]-=fx;
		(F_y)[a]-=fy;
		(F_x)[b]-=fx;
		(F_y)[b]-=fy;
	}
}

void BertaultLayout::r_Calc_On_Edge(node *v, edge *e, GraphAttributes &AG)
{
	node a=(*e)->source();
	node b=(*e)->target();
	int s=0;
	double x_diff=i.x-AG.x(*v);
	double y_diff=i.y-AG.y(*v);

	//determines the section in which the line-segment (v,i) lies
	if(x_diff>=0)
	{
		if(y_diff>=0)
		{
			if(x_diff>=y_diff)
				s=1;
			else
				s=2;
		}
		else
		{
			if(x_diff>=-(y_diff))
				s=8;
			else
				s=7;
		}
	}
	else
	{
		if(y_diff>=0)
		{
			if(-(x_diff)>=y_diff)
				s=4;
			else
				s=3;
		}
		else
		{
			if(-(x_diff)>=-(y_diff))
				s=5;
			else
				s=6;
		}
	}

	OGDF_ASSERT(s!=0);			//section>=1
	double max_radius=(sqrt(x_diff*x_diff+y_diff*y_diff))/3;

	//cout << "node:" << (*v)->index() << "edge:" << (*e)->index() << "between" << a->index() << "and" << b->index() << "INSIDE\nsection:" << s << "x_a:" << AG.x(a) << "x_b:" << AG.x(b) << "y_a:" << AG.y(a) << "y_b:" << AG.y(b) << "i_x:" << i.x << "i_y:" << i.y << "dist_v-i:" << max_radius*3 << "\n\n" ;

  //  if(max_radius>=0)
	//{
		int r,num;
		//determines which sections should have their values changed
	for(r=s-2;r<=(s+2);r++)
	{
			num=1+((r-1)%8);
			if(num<=0)
				num+=8;
			(sect)[*v].R[num]=min((sect)[*v].R[num],max_radius);
		}
		for(r=s+2;r<=(s+6);r++)
	{
			 num=1+((r-1)%8);
			 if(num<=0)
				 num+=8;
			(sect)[a].R[num]=min((sect)[a].R[num],max_radius);
			(sect)[b].R[num]=min((sect)[b].R[num],max_radius);
		}
	//}
	/*int i;
	for(i=1;i<9;i++)
	{
		cout << sect[*v].R[i];
		cout << " ";
	}
	cout << "\n";*/
}


void BertaultLayout::r_Calc_Outside_Edge(node *v, edge *e, GraphAttributes &AG)
{
	node a=(*e)->source();
	node b=(*e)->target();
	double dav=sqrt((AG.x(*v)-AG.x(a))*(AG.x(*v)-AG.x(a))+(AG.y(*v)-AG.y(a))*(AG.y(*v)-AG.y(a)));
	double dbv=sqrt((AG.x(*v)-AG.x(b))*(AG.x(*v)-AG.x(b))+(AG.y(*v)-AG.y(b))*(AG.y(*v)-AG.y(b)));
	//cout << "node:" << (*v)->index() << "edge:" << (*e)->index() << "between" << a->index() << "and" << b->index() << "OUTSIDE\ndist_v_a:" << dav << "dist_v_b:" << dbv << "\n\n" ;

	int r;
	for(r=1;r<=8;r++)
	{
		(sect)[*v].R[r]=min((sect)[*v].R[r],min(dav,dbv)/3);
		(sect)[a].R[r]=min((sect)[a].R[r],dav/3);
		(sect)[b].R[r]=min((sect)[b].R[r],dbv/3);
	}

	/*int i;
	for(i=1;i<9;i++)
	{
		cout << sect[*v].R[i];
		cout << " ";
	}
	cout << "\n";*/
}

void BertaultLayout::move(node *v, GraphAttributes &AG)
{
	int s=0;
	double x_diff=(F_x)[*v];
	double y_diff=(F_y)[*v];

	//determines the section in which the node has to move
	if(x_diff>=0)
	{
		if(y_diff>=0)
		{
			if(x_diff>=y_diff)
				s=1;
			else
				s=2;
		}
		else
		{
			if(x_diff>=-(y_diff))
				s=8;
			else
				s=7;
		}
	}
	else
	{
		if(y_diff>=0)
		{
			if(-(x_diff)>=y_diff)
				s=4;
			else
				s=3;
		}
		else
		{
			if(-(x_diff)>=-(y_diff))
				s=5;
			else
				s=6;
		}
	}

	OGDF_ASSERT(s!=0);

	double mov_mag=sqrt((x_diff*x_diff)+(y_diff*y_diff));		// the length of the move
	if((sect)[*v].R[s]<mov_mag)									// if move is greater than zone(section) radius
	{															// magnitudes of forces are normalised so that the length becomes equal to the radius
		(F_x)[*v]=((F_x)[*v]/mov_mag)*(sect)[*v].R[s];					// and the move will now take the node on the arc of that section
		(F_y)[*v]=((F_y)[*v]/mov_mag)*(sect)[*v].R[s];
	}
	//moves the node
	AG.x(*v)+=(F_x)[*v];
	AG.y(*v)+=(F_y)[*v];
}

void BertaultLayout::initPositions(GraphAttributes &AG, char c)
{
	if( (AG.attributes() & GraphAttributes::nodeGraphics) == 0 && (c=='c'||c=='m'||c=='r'))
	{
		if(req_length==0)
			req_length=50;
		AG.initAttributes((AG.attributes()|GraphAttributes::nodeGraphics|GraphAttributes::edgeGraphics|GraphAttributes::nodeStyle|GraphAttributes::edgeStyle));
		const Graph &G = AG.constGraph();
		int m = (int) sqrt((double)G.numberOfNodes());
		int cnth=0,cntc=0;
		int dim = (int)(req_length*G.numberOfNodes()/2);
		node v;
		srand ( (unsigned int) time(NULL) );
		forall_nodes(v,G)
		{
			if(c=='r')
			{
				int flag=1;
				while(flag==1)
				{
					AG.x(v)=(double)(rand()%dim)-dim/2;
					AG.y(v)=(double)(rand()%dim)-dim/2;
					flag=0;
					node x;
					forall_nodes(x,G)
					{
						if(x==v)
							break;
						if(AG.x(v)==AG.x(x)&&AG.y(v)==AG.y(x))
						{
							flag=1;
							break;
						}
					}
				}
			}
			else
			{
				int flag=1;
				while(flag==1)
				{
					if(c=='c')
					{
						double r=req_length*(cntc+1)/2;
						double ang=(2*Math::pi/m)*cnth;
						double cs=cos(ang);
						double sn=sin(ang);
						if((cs<1.0e-8&&cs>0)||(cs>-1.0e-8&&cs<0))
						{
							if(sn<0)
								sn=-1;
							else
								sn=1;

							cs=0;
						}
						if((sn<1.0e-8&&sn>0)||(sn>-1.0e-8&&sn<0))
						{
							if(cs<0)
								cs=-1;
							else
								cs=1;

							sn=0;
						}

						AG.x(v)=r*cs;
						AG.y(v)=r*sn;
						/*  if((AG.x(v)*AG.x(v)+AG.y(v)*AG.y(v))!=r*r)
						{
							cout << (AG.x(v)*AG.x(v)+AG.y(v)*AG.y(v))-(r*r);
							cout << "\n";
						}*/
					}
					else if(c=='m')
					{
						AG.x(v)=req_length*cnth/2;
						AG.y(v)=req_length*cntc/2;
					}

					node x;
					flag=0;
					forall_nodes(x,G)
					{
						if(x==v)
							break;
						if(AG.x(v)==AG.x(x)&&AG.y(v)==AG.y(x))
						{
							flag=1;
							cnth--;
							break;
						}
					}
				}

				cnth++;
				if(cnth==m)
				{
					cnth=0;
					cntc++;
				}
			}

			AG.width(v)=req_length/10;
			AG.height(v)=req_length/10;
		}
	}
}

void BertaultLayout::preprocess(GraphAttributes &AG)
{
	node v;
	const Graph &G = AG.constGraph();
	surr.init(0,G.numberOfNodes()-1,0,G.numberOfEdges()-1,0);
	GraphCopy G1(G);
	GraphAttributes AG1(G1);
	AG1.setDirected(AG.directed());
	forall_nodes(v,G1)
	{
		AG1.x(v)=AG.x(G1.original(v));
		AG1.y(v)=AG.y(G1.original(v));
		AG1.width(v)=AG.width(G1.original(v));
		AG1.height(v)=AG.height(G1.original(v));
	}

	labelling(AG1);
	crossingPlanarize(AG1);


	/*surr.init(G);
	forall_nodes(v,G)
	{
		surr[v].init(G);
		forall_edges(e,G)
		{
			surr[v][e]=0;
		}
	}*/
//    getOriginal(AG);

	PlanRep PG(AG1);
	int numCC = PG.numberOfCCs();
	//GraphAttributes PAG(PG, GraphAttributes::nodeGraphics|GraphAttributes::edgeGraphics);
	int i;
	//List< node > list;

	//cout << numCC <<"\n";

	List<CCElement*> forest;
	Array<CCElement> Carr(numCC);

	for(i = 0; i < numCC; i++)
	{
		Carr[i].init(i);
		Carr[i].faceNum=-1;
	}

	for(i = 0; i < numCC; i++)
	{
		CCElement* new1=&Carr[i];
		int rootnum=0,flag=0;
			while(rootnum<forest.size())
			{
				/* //CCElement *Croot=&(**(forest.get(rootnum)));
				cout << "Children of " << (**(forest.get(rootnum))).num <<"bfor calling :\n";
				int k;
				for(k=0;k<(**(forest.get(rootnum))).child.size();k++)
				{
					cout << (*((**(forest.get(rootnum))).child.get(k)))->num << "\n";
				}
				cout << "Inserting " << (*new1).num << " into " << (**(forest.get(rootnum))).num << "\n";
				//*/

				int retv=insert(new1,&(**(forest.get(rootnum))),AG1,PG);
				if(retv==2)
				{
					flag=1;
					break;
				}
				else if(retv==1)
				{
					(**(forest.get(rootnum))).root=false;
					//ListIterator<CCElement> l((*(forest.get(rootnum))));
					//cout << (**forest.get(rootnum)).num << " deleted\n";
					forest.del(forest.get(rootnum));
					rootnum--;
				}
				rootnum++;
			}

		if(flag==0)
		{
			(*new1).faceNum=-1;
			(*new1).root=true;

			forest.pushBack(&(*new1));
			/*
			cout << (*new1).num << " is a root now " << "\nChildren of " << (*new1).num << ":\n";
			int j;
			for(j=0;j<(*new1).child.size();j++)
			{
				cout << (*((*new1).child.get(j)))->num << "\n";
			}
			//*/
		}
	}

	// Uncomment below statements to see output... for debugging use


	//cout << "forest size : " << forest.size();
/*
	for(i=0;i<forest.size();i++)
	{
		cout << "\nRoot : " << i <<"\n";
		compute(*forest.get(i),PG,AG1,G1);
	}
//*/

/*
	forall_nodes(v,G)
	{
		forall_edges(e,G)
		{
			cout<< surr(v->index(),e->index()) <<" ";
		}
		cout << "\n";
	}
	//*/
/*
	node n=G.chooseNode();
	AG.fillColor(n)="RED";
	forall_edges(e,G)
	{
		if(surr(n->index(),e->index())==true)
			AG.strokeColor(e)="RED";
	}

	GraphIO::writeGML(AG, "planarized.gml");
	//*/
}

void BertaultLayout::labelling(GraphAttributes &AG)
{
	AG.initAttributes(GraphAttributes::edgeIntWeight);
	edge e;
	forall_edges(e,AG.constGraph())
	{
		AG.intWeight(e)=e->index();
	}
}

void BertaultLayout::crossingPlanarize(GraphAttributes &AG)
{
	Graph &G= const_cast<Graph&> (AG.constGraph());

	edge e;
	forall_edges(e,G)
	{
		edge i;
		forall_rev_edges(i,G)
		{
			if(i==e)
				break;

			node a=e->source();
			node b=e->target();
			double m=(AG.y(a)-AG.y(b))/(AG.x(a)-AG.x(b));
			double c=AG.y(a)-m*AG.x(a);

			node x=i->source();
			node y=i->target();

			if(a!=x&&a!=y&&b!=x&&b!=y)
			{
				double m2=(AG.y(x)-AG.y(y))/(AG.x(x)-AG.x(y));
				double c2=AG.y(x)-m2*AG.x(x);

				double ainc=(AG.y(a)-m2*AG.x(a)-c2),binc=(AG.y(b)-m2*AG.x(b)-c2),xinc=AG.y(x)-m*AG.x(x)-c,yinc=AG.y(y)-m*AG.x(y)-c;
				//int temp=AG.intWeight(e);

				if(xinc*yinc<0&&ainc*binc<0)
				{
					//cout << "edge " << e->index() << " and edge " << i->index() << " between nodes " << a->index() << " and " << b->index() << " and nodes " << x->index() << " and " << y->index() << "\n";
					int temp=AG.intWeight(e);
					edge enew=G.split(e);
					node nnew=enew->source();
					AG.width(nnew)=AG.width(a);
					AG.height(nnew)=AG.height(a);
					AG.x(nnew)=(c2-c)/(m-m2);
					AG.y(nnew)=m*AG.x(nnew)+c;
					AG.intWeight(enew)=temp;
					edge xn=G.newEdge(x,nnew);
					AG.intWeight(xn)=AG.intWeight(i);
					AG.intWeight(G.newEdge(nnew,y))=AG.intWeight(i);
					G.delEdge(i);
				}
			}
			//*/
		}
	}
}

int BertaultLayout::insert(CCElement *new1,CCElement *element,GraphAttributes &PAG,PlanRep &PG)
{
	int contface=contained(new1,element,PAG,PG);
	if(contface!=-1)
	{
		int flag=0,i;
	/*
		cout << "Children of " << (*element).num <<"bfor :\n";
		for(i=0;i<(*element).child.size();i++)
		{
			cout << (*((*element).child.get(i)))->num << "\n";
		}
	//*/
		if((*element).child.size()!=0)
		{

			for(i=0;i<(*element).child.size();i++)
			{
				CCElement *child=&(**((*element).child.get(i)));
				if(child->faceNum==contface)
				{
					//cout << "Gonna insert " << (*new1).num << " in " << child->num << "\n";

					int retv=insert(new1,child,PAG,PG);
					if(retv==2)
					{
						flag=1;
						break;
					}
					else if(retv==1)
					{
						i--;
					}
				}
			}
		}

		if(flag==0)
		{
			(*new1).parent=&(*element);
			(*new1).faceNum=contface;
			(*element).child.pushBack(&(*new1));
	  /*
			cout << (*new1).num << " is child of " << new1->parent->num << "\nChildren of " << (*element).num << ":\n";
			for(i=0;i<(*element).child.size();i++)
			{
				cout << (*((*element).child.get(i)))->num << "\n";
			}

			cout << "Children of " << (*new1).num <<" ("<< (*new1).child.size() <<"):\n";
			for(i=0;i<(*new1).child.size();i++)
			{
				cout << (*((*new1).child.get(i)))->num << "\n";
			}
		//*/

		}
		return 2;
	}
	else
	{
		contface=contained(element,new1,PAG,PG);
		if(contface!=-1)
		{
			//ListIterator<CCElement*> l;
			if(!(*element).root)
			{
				//cout << "deleting " << (**(*element).parent->child.get((*element).parent->child.search(&(*element)))).num << "\n";
				(*element).parent->child.del((*element).parent->child.get((*element).parent->child.search(&(*element))));
			}
		   /*
			else
				cout << (*element).num << " is a root :(\n";
			//*/
			(*element).faceNum=contface;
			(*element).parent=&(*new1);
			(*new1).child.pushBack(&(*element));

			/*
			cout << "Pushed " << (*element).num << " into " << (*new1).num << "\n";

			cout << (*new1).num << " is parent of " << (*element).num << "\nChildren of " << (*new1).num << ":\n";

			int i;
			for(i=0;i<(*new1).child.size();i++)
			{
				cout << (**((*new1).child.get(i))).num << "\n";
			}
			//*/
			return 1;
		}
		else
		{
			return 0;
		}
	}
}


int BertaultLayout::contained(CCElement* new1,CCElement* element,GraphAttributes &PAG,PlanRep &PG)
{
	//PlanRep &PG= (PlanRep&)(const_cast<Graph&> (PAG.constGraph()));

	PG.initCC(new1->num);
	node v;

	v=PG.chooseNode();
	double yc=PAG.y(PG.original(v));
	double xc=PAG.x(PG.original(v));

	//cout << "Number of edges in " << new1->num << " is " << PG.numberOfEdges() << "\n";

	PG.initCC(element->num);
	//cout << "Number of edges in " << element->num << " is " << PG.numberOfEdges() << "\n";
	//cout << "Chosen " << xc;
	ConstCombinatorialEmbedding E(PG);
	E.computeFaces();
	face f;
	forall_faces(f,E)
	{
		int crossings=0;
		adjEntry adj;
		List< int > edges;
		forall_face_adj(adj,f)
		{

			if(edges.search(adj->theEdge()->index())==-1)
			{
				edges.pushBack(adj->theEdge()->index());
				node x=adj->theEdge()->source();
				node y=adj->theEdge()->target();
				double m= (PAG.y(PG.original(x))-PAG.y(PG.original(y)))/(PAG.x(PG.original(x))-PAG.x(PG.original(y)));

				double c= PAG.y(PG.original(x))-m*PAG.x(PG.original(x));
				if((PAG.y(PG.original(x))-yc)*(PAG.y(PG.original(y))-yc)<=0&&((yc-c)/m)>=xc)
				{
					crossings++;
				}

			}

		}

		if(crossings%2!=0)
		{
			//cout << "yo";
			return f->index();
		}
	}
	return -1;
}


void BertaultLayout::compute(CCElement* element,PlanRep &PG,GraphAttributes &AG1,GraphCopy &G1)
{
	int num=element->num;
	PG.initCC(num);
	ConstCombinatorialEmbedding E(PG);
	//E.computeFaces();
	//FaceArray< List<edge> > farray(E);
	face f;
	forall_faces(f,E)
	{
		//cout<<"yeah";
		adjEntry adj;
		forall_face_adj(adj,f)
		{
			node ver=adj->theNode();
			node ver2=adj->twinNode();
			bool dum=false,dum2=false;
			if(G1.isDummy(PG.original(ver)))
				dum=true;
			if(G1.isDummy(PG.original(ver2)))
				dum2=true;

			node v=G1.original(PG.original(ver));
			node v2=G1.original(PG.original(ver2));

			//cout << v->index() << "\n";
			adjEntry adj1 = f->firstAdj(), adj3 = adj1;
			do {
			//cout<<"lala\n";
			if(!dum)
				surr(v->index(),AG1.intWeight(PG.original(adj3->theEdge())))=true;
			if(!dum2)
				surr(v2->index(),AG1.intWeight(PG.original(adj3->theEdge())))=true;
			adj3 = adj3->faceCycleSucc();
			} while (adj3->index() != adj1->index());

			int j;
			for(j=0;j<element->child.size();j++)
			{
				if((*(element->child.get(j)))->faceNum==f->index())
				{
					PG.initCC((*(element->child.get(j)))->num);
					ConstCombinatorialEmbedding E2(PG);
					E2.computeFaces();
					face f2;
					adjEntry adj2;
					forall_faces(f2,E2)
					{
						forall_face_adj(adj2,f2)
						{
							if(!dum)
								surr(v->index(),AG1.intWeight(PG.original(adj2->theEdge())))=true;
							if(!dum2)
								surr(v2->index(),AG1.intWeight(PG.original(adj2->theEdge())))=true;
						}
					}
				}
			}

			if(element->faceNum!=-1)
			{
				PG.initCC(element->parent->num);
				ConstCombinatorialEmbedding E3(PG);
				E3.computeFaces();
				face f2;
				//int flag=0;
				forall_faces(f2,E3)
				{
					if(f2->index()==element->faceNum)
					{
						//flag=1;
						break;
					}
				}

				adjEntry adj2;
				forall_face_adj(adj2,f2)
				{
					if(!dum)
						surr(v->index(),AG1.intWeight(PG.original(adj2->theEdge())))=true;
					if(!dum2)
						surr(v2->index(),AG1.intWeight(PG.original(adj2->theEdge())))=true;
				}
			}
			PG.initCC(num);
		}
	}
	//cout << num <<" Done\n";

	int i;
	for(i=0;i<element->child.size();i++)
	{
		compute(*(element->child.get(i)),PG,AG1,G1);
	}

}

int BertaultLayout::edgeCrossings(GraphAttributes &AG)
{
	const Graph &G = AG.constGraph();
	int crossings=0;
	edge e;
	forall_edges(e,G)
	{
		node a=e->source();
		node b=e->target();
		double m=(AG.y(a)-AG.y(b))/(AG.x(a)-AG.x(b));
		double c=AG.y(a)-m*AG.x(a);
		edge i;
		forall_rev_edges(i,G)
		{
			if(i==e)
				break;

			node x=i->source();
			node y=i->target();
			double m2=(AG.y(x)-AG.y(y))/(AG.x(x)-AG.x(y));
			double c2=AG.y(x)-m2*AG.x(x);
			double distab=sqrt((AG.x(a)-AG.x(b))*(AG.x(a)-AG.x(b))+(AG.y(a)-AG.y(b))*(AG.y(a)-AG.y(b)));
			double distxy=sqrt((AG.x(x)-AG.x(y))*(AG.x(x)-AG.x(y))+(AG.y(x)-AG.y(y))*(AG.y(x)-AG.y(y)));
			double distax=sqrt((AG.x(a)-AG.x(x))*(AG.x(a)-AG.x(x))+(AG.y(a)-AG.y(x))*(AG.y(a)-AG.y(x)));
			double distay=sqrt((AG.x(a)-AG.x(y))*(AG.x(a)-AG.x(y))+(AG.y(a)-AG.y(y))*(AG.y(a)-AG.y(y)));
			double distbx=sqrt((AG.x(x)-AG.x(b))*(AG.x(x)-AG.x(b))+(AG.y(x)-AG.y(b))*(AG.y(x)-AG.y(b)));
			double distby=sqrt((AG.x(y)-AG.x(b))*(AG.x(y)-AG.x(b))+(AG.y(y)-AG.y(b))*(AG.y(y)-AG.y(b)));
			double d=distab+distxy;

			if(a!=x&&a!=y&&b!=x&&b!=y)
			{
				double ainc=(AG.y(a)-m2*AG.x(a)-c2),binc=(AG.y(b)-m2*AG.x(b)-c2),xinc=AG.y(x)-m*AG.x(x)-c,yinc=AG.y(y)-m*AG.x(y)-c;

				if(((xinc*yinc<0&&ainc*binc<0)||(xinc*yinc==0&&ainc*binc<0)||(xinc*yinc<0&&ainc*binc==0)))
				{
				 // cout << "edge " << e->index() << " and edge " << i->index() << " between nodes " << a->index() << " and " << b->index() << " and nodes " << x->index() << " and " << y->index() << "\n";
					crossings++;
				}
				else
				{
					if(m==m2&&c==c2&&distax<d&&distay<d&&distbx<d&&distby<d)
					{
						//cout << "edge " << e->index() << " and edge " << i->index() << " between nodes " << a->index() << " and " << b->index() << " and nodes " << x->index() << " and " << y->index() << "OVERLAP\n";
						crossings+=2;
					}
				}
			}
			else if(m==m2&&c==c2&&distax<d&&distay<d&&distbx<d&&distby<d&&((a!=y&&b!=x&&b!=y)||(a!=x&&b!=x&&b!=y)||(a!=y&&a!=x&&b!=y)||(a!=y&&a!=x&&b!=x)))
			{
				//cout << "edge " << e->index() << " and edge " << i->index() << " between nodes " << a->index() << " and " << b->index() << " and nodes " << x->index() << " and " << y->index() << "OVERLAP\n";
				crossings+=1;
			}

		}
	}
	return crossings;
}

double BertaultLayout::edgelength(GraphAttributes &AG)
{
	EdgeArray<double> el;
	const Graph &G = AG.constGraph();
	el.init(G);
	edge e;
	double mean=0,stdev=0;
	forall_edges(e,G)
	{
		node a=e->source();
		node b=e->target();
		el[e]=sqrt((AG.x(a)-AG.x(b))*(AG.x(a)-AG.x(b))+(AG.y(a)-AG.y(b))*(AG.y(a)-AG.y(b)));
		mean+=el[e];
	}
	mean=mean/(G.numberOfEdges());
	forall_edges(e,G)
	{
		stdev+=(el[e]-mean)*(el[e]-mean);
	}
	stdev=sqrt(stdev/(G.numberOfEdges()))/mean;
   /* cout << "\n";
	cout << "mean : ";
	cout << mean;
	cout << "\n";*/
	return stdev;
}

double BertaultLayout::nodeDistribution(GraphAttributes &AG)
{
	const Graph &G = AG.constGraph();
	if(G.numberOfNodes()<2)
		return -1;
	node v = G.firstNode();
	double
	  minx = AG.x(v),
	  maxx = AG.x(v),
	  miny = AG.y(v),
	  maxy = AG.y(v);
	for (v = v->succ(); v; v = v->succ()) {
		if (AG.x(v) > maxx)
			maxx = AG.x(v);
		if (AG.x(v) < minx)
			minx = AG.x(v);
		if (AG.y(v) > maxy)
			maxy = AG.y(v);
		if (AG.y(v) < miny)
			miny = AG.y(v);
	}

	int rows=8,columns=8,i,j;
	double sizex = (maxx-minx)/(double)(columns-1), sizey = (maxy-miny)/(double)(columns-1);
	double startx = minx-sizex/2, /* endx = maxx+sizex/2, */ starty = miny-sizey/2 /* , endy = maxy+sizey/2 */ ;
	//int box[rows][columns];
	Array2D<int> box(0,rows-1,0,columns-1);

	for(i=0; i<rows; i++)
		for(j=0; j<columns; j++)
			box(i,j) = 0;
	if(maxy != miny && maxx != minx)
	{
		forall_nodes(v,G)
		{
			box((int)((AG.y(v)-starty)/sizey), (int)((AG.x(v)-startx)/sizex))++;
		}

		double mean=(double)G.numberOfNodes()/(double)(rows*columns);
		double stdev = 0;
		for(i=0;i<rows;i++)
			for(j=0;j<columns;j++)
			{
			  //  cout << box(i,j);
				stdev+=((double)(box(i,j))-mean)*((double)(box(i,j))-mean);
			}

	   // cout << "\n";
		stdev=sqrt(stdev/(rows*columns))/mean;
		return stdev;
	}
	else
		return -1;
}

}//namespace ogdf
