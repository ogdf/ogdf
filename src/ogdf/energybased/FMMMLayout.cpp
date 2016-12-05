/** \file
 * \brief Implementation of Fast Multipole Multilevel Method (FM^3).
 *
 * \author Stefan Hachul
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.md in the OGDF root directory for details.
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
 * License along with this program; if not, see
 * http://www.gnu.org/copyleft/gpl.html
 */


#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/basic/Math.h>
#include <ogdf/internal/energybased/numexcept.h>
#include <ogdf/internal/energybased/MAARPacking.h>
#include <ogdf/internal/energybased/Multilevel.h>
#include <ogdf/internal/energybased/Rectangle.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/basic.h>

#include <cmath>

namespace ogdf {


FMMMLayout::FMMMLayout()
{
	initialize_all_options();
}


//--------------------------- most important functions --------------------------------

void FMMMLayout::call(GraphAttributes &GA)
{
	const Graph &G = GA.constGraph();
	EdgeArray<double> edgelength(G,1.0);
	call(GA,edgelength);
}

void FMMMLayout::call(ClusterGraphAttributes &GA)
{
	const Graph &G = GA.constGraph();
	//compute depth of cluster tree, also sets cluster depth values
	const ClusterGraph &CG = GA.constClusterGraph();
	int cdepth = CG.treeDepth();
	EdgeArray<double> edgeLength(G);
	//compute lca of end vertices for each edge
	for(edge e : G.edges)
	{
		edgeLength[e] = cdepth - CG.clusterDepth(CG.commonCluster(e->source(),e->target())) + 1;
		OGDF_ASSERT(edgeLength[e] > 0);
	}
	call(GA,edgeLength);
	GA.updateClusterPositions();
}


void FMMMLayout::call(GraphAttributes &GA, const EdgeArray<double> &edgeLength)
{
	const Graph &G = GA.constGraph();
	NodeArray<NodeAttributes> A(G);       //stores the attributes of the nodes (given by L)
	EdgeArray<EdgeAttributes> E(G);       //stores the edge attributes of G
	Graph G_reduced;                      //stores a undirected simple and loopfree copy
										//of G
	EdgeArray<EdgeAttributes> E_reduced;  //stores the edge attributes of G_reduced
	NodeArray<NodeAttributes> A_reduced;  //stores the node attributes of G_reduced

	if(G.numberOfNodes() > 1)
	{
		GA.clearAllBends();//all edges are straight-line
		if(useHighLevelOptions())
			update_low_level_options_due_to_high_level_options_settings();
		import_NodeAttributes(G,GA,A);
		import_EdgeAttributes(G,edgeLength,E);

		double t_total;
		usedTime(t_total);
		max_integer_position = pow(2.0,maxIntPosExponent());
		init_ind_ideal_edgelength(G,A,E);
		make_simple_loopfree(G,A,E,G_reduced,A_reduced,E_reduced);
		call_DIVIDE_ET_IMPERA_step(G_reduced,A_reduced,E_reduced);
		if(allowedPositions() != apAll)
			make_positions_integer(G_reduced,A_reduced);
		time_total = usedTime(t_total);

		export_NodeAttributes(G_reduced,A_reduced,GA);
	}
	else //trivial cases
	{
		if(G.numberOfNodes() == 1 )
		{
			node v = G.firstNode();
			GA.x(v) = 0;
			GA.y(v) = 0;
		}
	}
}


void FMMMLayout::call(GraphAttributes &AG, char* ps_file)
{
	call(AG);
	create_postscript_drawing(AG,ps_file);
}


void FMMMLayout::call(
	GraphAttributes &AG,
	const EdgeArray<double> &edgeLength,
	char* ps_file)
{
	call(AG,edgeLength);
	create_postscript_drawing(AG,ps_file);
}


void FMMMLayout::call_DIVIDE_ET_IMPERA_step(
	Graph& G,
	NodeArray<NodeAttributes>& A,
	EdgeArray<EdgeAttributes>& E)
{
	NodeArray<int> component(G); //holds for each node the index of its component
	number_of_components = connectedComponents(G,component);//calculate components of G
	Graph* G_sub = new Graph[number_of_components];
	NodeArray<NodeAttributes>* A_sub = new NodeArray<NodeAttributes>[number_of_components];
	EdgeArray<EdgeAttributes>* E_sub = new EdgeArray<EdgeAttributes>[number_of_components];
	create_maximum_connected_subGraphs(G,A,E,G_sub,A_sub,E_sub,component);

	if(number_of_components == 1)
		call_MULTILEVEL_step_for_subGraph(G_sub[0],A_sub[0],E_sub[0],-1);
	else
		for(int i = 0; i < number_of_components;i++)
			call_MULTILEVEL_step_for_subGraph(G_sub[i],A_sub[i],E_sub[i],i);

	pack_subGraph_drawings (A,G_sub,A_sub);
	delete_all_subGraphs(G_sub,A_sub,E_sub);
}


void FMMMLayout::call_MULTILEVEL_step_for_subGraph(
	Graph& G,
	NodeArray<NodeAttributes>& A,
	EdgeArray<EdgeAttributes>& E,
	int comp_index)
{
	energybased::Multilevel Mult;

	int max_level = 30;//sufficient for all graphs with upto pow(2,30) nodes!
	//adapt mingraphsize such that no levels are created beyond input graph.
	if (m_singleLevel) m_minGraphSize = G.numberOfNodes();
	Array<Graph*> G_mult_ptr(max_level+1);
	Array<NodeArray<NodeAttributes>*> A_mult_ptr (max_level+1);
	Array<EdgeArray<EdgeAttributes>*> E_mult_ptr (max_level+1);

	Mult.create_multilevel_representations(G,A,E,randSeed(),
				galaxyChoice(),minGraphSize(),
				randomTries(),G_mult_ptr,A_mult_ptr,
				E_mult_ptr,max_level);

	for(int i = max_level;i >= 0;i--)
	{
		if(i == max_level)
			create_initial_placement(*G_mult_ptr[i],*A_mult_ptr[i]);
		else
		{
			Mult.find_initial_placement_for_level(i,initialPlacementMult(),G_mult_ptr, A_mult_ptr,E_mult_ptr);
			update_boxlength_and_cornercoordinate(*G_mult_ptr[i],*A_mult_ptr[i]);
		}
		call_FORCE_CALCULATION_step(*G_mult_ptr[i],*A_mult_ptr[i],*E_mult_ptr[i],
				 i,max_level);
	}
	Mult.delete_multilevel_representations(G_mult_ptr,A_mult_ptr,E_mult_ptr,max_level);
}


void FMMMLayout::call_FORCE_CALCULATION_step(
	Graph& G,
	NodeArray<NodeAttributes>&A,
	EdgeArray<EdgeAttributes>& E,
	int act_level,
	int max_level)
{
	const int ITERBOUND = 10000;//needed to guarantee termination if
							 //stopCriterion() == scThreshold
	if(G.numberOfNodes() > 1)
	{
		int iter = 1;
		int max_mult_iter = get_max_mult_iter(act_level,max_level,G.numberOfNodes());
		double actforcevectorlength = threshold() + 1;

		NodeArray<DPoint> F_rep(G); //stores rep. forces
		NodeArray<DPoint> F_attr(G); //stores attr. forces
		NodeArray<DPoint> F (G); //stores resulting forces
		NodeArray<DPoint> last_node_movement(G);//stores the force vectors F of the last
												//iterations (needed to avoid oscillations)

		set_average_ideal_edgelength(G,E);//needed for easy scaling of the forces
		make_initialisations_for_rep_calc_classes(G);

		while( ((stopCriterion() == scFixedIterations)&&(iter <= max_mult_iter)) ||
			((stopCriterion() == scThreshold)&&(actforcevectorlength >= threshold())&&
			(iter <= ITERBOUND)) ||
			((stopCriterion() == scFixedIterationsOrThreshold)&&(iter <= max_mult_iter) &&
			(actforcevectorlength >= threshold())) )
		{//while
			calculate_forces(G,A,E,F,F_attr,F_rep,last_node_movement,iter,0);
			if(stopCriterion() != scFixedIterations)
				actforcevectorlength = get_average_forcevector_length(G,F);
			iter++;
		}//while

		if(act_level == 0)
			call_POSTPROCESSING_step(G,A,E,F,F_attr,F_rep,last_node_movement);

		deallocate_memory_for_rep_calc_classes();
	}
}


void FMMMLayout::call_POSTPROCESSING_step(
	Graph& G,
	NodeArray<NodeAttributes>& A,
	EdgeArray<EdgeAttributes>& E,
	NodeArray<DPoint>& F,
	NodeArray<DPoint>& F_attr,
	NodeArray<DPoint>& F_rep,
	NodeArray<DPoint>& last_node_movement)
{
	for(int i = 1; i<= 10; i++)
		calculate_forces(G,A,E,F,F_attr,F_rep,last_node_movement,i,1);

	if(resizeDrawing())
	{
		adapt_drawing_to_ideal_average_edgelength(G,A,E);
		update_boxlength_and_cornercoordinate(G,A);
	}

	for(int i = 1; i<= fineTuningIterations(); i++)
		calculate_forces(G,A,E,F,F_attr,F_rep,last_node_movement,i,2);

	if(resizeDrawing())
		adapt_drawing_to_ideal_average_edgelength(G,A,E);
}


//------------------------- functions for pre/post-processing -------------------------

void FMMMLayout::initialize_all_options()
{
	//setting high level options
	useHighLevelOptions(false);
	pageFormat(pfSquare);
	unitEdgeLength(LayoutStandards::defaultNodeSeparation());
	newInitialPlacement(false);
	qualityVersusSpeed(qvsBeautifulAndFast);

	//setting low level options
	//setting general options
	randSeed(100);
	edgeLengthMeasurement(elmBoundingCircle);
	allowedPositions(apInteger);
	maxIntPosExponent(40);

	//setting options for the divide et impera step
	pageRatio(1.0);
	stepsForRotatingComponents(10);
	tipOverCCs(toNoGrowingRow);
	minDistCC(LayoutStandards::defaultCCSeparation());
	presortCCs(psDecreasingHeight);

	//setting options for the multilevel step
	minGraphSize(50);
	galaxyChoice(gcNonUniformProbLowerMass);
	randomTries(20);
	maxIterChange(micLinearlyDecreasing);
	maxIterFactor(10);
	initialPlacementMult(ipmAdvanced);
	m_singleLevel = false;

	//setting options for the force calculation step
	forceModel(fmNew);
	springStrength(1);
	repForcesStrength(1);
	repulsiveForcesCalculation(rfcNMM);
	stopCriterion(scFixedIterationsOrThreshold);
	threshold(0.01);
	fixedIterations(30);
	forceScalingFactor(0.05);
	coolTemperature(false);
	coolValue(0.99);
	initialPlacementForces(ipfRandomRandIterNr);

	//setting options for postprocessing
	resizeDrawing(true);
	resizingScalar(1);
	fineTuningIterations(20);
	fineTuneScalar(0.2);
	adjustPostRepStrengthDynamically(true);
	postSpringStrength(2.0);
	postStrengthOfRepForces(0.01);

	//setting options for different repulsive force calculation methods
	frGridQuotient(2);
	nmTreeConstruction(rtcSubtreeBySubtree);
	nmSmallCell(scfIteratively);
	nmParticlesInLeaves(25);
	nmPrecision(4);
}


void FMMMLayout::update_low_level_options_due_to_high_level_options_settings()
{
	PageFormatType pf = pageFormat();
	double uel = unitEdgeLength();
	bool nip = newInitialPlacement();
	QualityVsSpeed qvs = qualityVersusSpeed();

	//update
	initialize_all_options();
	useHighLevelOptions(true);
	pageFormat(pf);
	unitEdgeLength(uel);
	newInitialPlacement(nip);
	qualityVersusSpeed(qvs);

	if(pageFormat() == pfSquare)
		pageRatio(1.0);
	else if(pageFormat() ==pfLandscape)
		pageRatio(1.4142);
	else //pageFormat() == pfPortrait
		pageRatio(0.7071);

	if(newInitialPlacement())
		initialPlacementForces(ipfRandomTime);
	else
		initialPlacementForces(ipfRandomRandIterNr);

	if(qualityVersusSpeed() == qvsGorgeousAndEfficient)
	{
		fixedIterations(60);
		fineTuningIterations(40);
		nmPrecision(6);
	}
	else if(qualityVersusSpeed() == qvsBeautifulAndFast)
	{
		fixedIterations(30);
		fineTuningIterations(20);
		nmPrecision(4);
	}
	else //qualityVersusSpeed() == qvsNiceAndIncredibleSpeed
	{
		fixedIterations(15);
		fineTuningIterations(10);
		nmPrecision(2);
	}
}


void FMMMLayout::import_NodeAttributes(
	const Graph& G,
	GraphAttributes& GA,
	NodeArray<NodeAttributes>& A)
{
	DPoint position;

	for(node v : G.nodes)
	{
		position.m_x = GA.x(v);
		position.m_y = GA.y(v);
		A[v].set_NodeAttributes(GA.width(v),GA.height(v),position,nullptr,nullptr);
	}
}


void FMMMLayout::import_EdgeAttributes(
	const Graph& G,
	const EdgeArray<double>& edgeLength,
	EdgeArray<EdgeAttributes>& E)
{
	double length;

	for(edge e : G.edges)
	{
		if(edgeLength[e] > 0) //no negative edgelength allowed
			length = edgeLength[e];
		else
			length = 1;

		E[e].set_EdgeAttributes(length,nullptr,nullptr);
	}
}


void FMMMLayout::init_ind_ideal_edgelength(
	const Graph& G,
	NodeArray<NodeAttributes>& A,
	EdgeArray<EdgeAttributes>& E)
{
	if (edgeLengthMeasurement() == elmMidpoint)
		for(edge e : G.edges)
			E[e].set_length(E[e].get_length() * unitEdgeLength());

	else //(edgeLengthMeasurement() == elmBoundingCircle)
	{
		set_radii(G,A);
		for(edge e : G.edges)
			E[e].set_length(E[e].get_length() * unitEdgeLength() + radius[e->source()]
				+ radius[e->target()]);
	}
}


void FMMMLayout::set_radii(const Graph& G, NodeArray<NodeAttributes>& A)
{
	radius.init(G);
	for(node v : G.nodes)
	{
		double w = A[v].get_width()/2;
		double h = A[v].get_height()/2;
		radius[v] = sqrt(w*w+ h*h);
	}
}


void FMMMLayout::export_NodeAttributes(
	Graph& G_reduced,
	NodeArray<NodeAttributes>& A_reduced,
	GraphAttributes& GA)
{
	for(node v_copy : G_reduced.nodes)
	{
		GA.x(A_reduced[v_copy].get_original_node()) =  A_reduced[v_copy].get_position().m_x;
		GA.y(A_reduced[v_copy].get_original_node()) =  A_reduced[v_copy].get_position().m_y;
	}
}


void FMMMLayout::make_simple_loopfree(
	const Graph& G,
	NodeArray<NodeAttributes>& A,
	EdgeArray<EdgeAttributes>E,
	Graph& G_reduced,
	NodeArray<NodeAttributes>& A_reduced,
	EdgeArray<EdgeAttributes>& E_reduced)
{
	//create the reduced Graph G_reduced and save in A/E links to node/edges of G_reduced
	//create G_reduced as a copy of G without selfloops!

	G_reduced.clear();
	for (node v_orig : G.nodes)
		A[v_orig].set_copy_node(G_reduced.newNode());

	for (edge e_orig : G.edges)
	{
		node u_orig = e_orig->source();
		node v_orig = e_orig->target();
		if (u_orig != v_orig)
			E[e_orig].set_copy_edge(G_reduced.newEdge(A[u_orig].get_copy_node(),
			A[v_orig].get_copy_node()));
		else
			E[e_orig].set_copy_edge(nullptr);//mark this edge as deleted
	}

	//remove parallel (and reversed) edges from G_reduced
	EdgeArray<double> new_edgelength(G_reduced);
	List<edge> S;
	S.clear();
	delete_parallel_edges(G, E, G_reduced, S, new_edgelength);

	//make A_reduced, E_reduced valid for G_reduced
	A_reduced.init(G_reduced);
	E_reduced.init(G_reduced);

	//import information for A_reduced, E_reduced and links to the original nodes/edges
	//of the copy nodes/edges
	for (node v_orig : G.nodes)
	{
		node v_reduced = A[v_orig].get_copy_node();
		A_reduced[v_reduced].set_NodeAttributes(A[v_orig].get_width(), A[v_orig].
			get_height(), A[v_orig].get_position(),
			v_orig, nullptr);
	}
	for (edge e_orig : G.edges)
	{
		edge e_reduced = E[e_orig].get_copy_edge();
		if (e_reduced != nullptr)
			E_reduced[e_reduced].set_EdgeAttributes(E[e_orig].get_length(), e_orig, nullptr);
	}

	//update edgelength of copy edges in G_reduced associated with a set of parallel
	//edges in G
	update_edgelength(S, new_edgelength, E_reduced);
}


void FMMMLayout::delete_parallel_edges(
	const Graph& G,
	EdgeArray<EdgeAttributes>& E,
	Graph& G_reduced,
	List<edge>& S,
	EdgeArray<double>& new_edgelength)
{
	energybased::EdgeMaxBucketFunc MaxSort;
	energybased::EdgeMinBucketFunc MinSort;
	energybased::Edge f_act;
	List<energybased::Edge> sorted_edges;
	EdgeArray<edge> original_edge (G_reduced); //helping array
	int save_s_index,save_t_index;
	int counter = 1;
	Graph* Graph_ptr = &G_reduced;

	//save the original edges for each edge in G_reduced
	for(edge e_act : G.edges)
		if(E[e_act].get_copy_edge() != nullptr) //e_act is no self_loops
			original_edge[E[e_act].get_copy_edge()] = e_act;

	for(edge e_act : G_reduced.edges)
	{
		f_act.set_Edge(e_act,Graph_ptr);
		sorted_edges.pushBack(f_act);
	}

	sorted_edges.bucketSort(0,G_reduced.numberOfNodes()-1,MaxSort);
	sorted_edges.bucketSort(0,G_reduced.numberOfNodes()-1,MinSort);

	edge e_save = nullptr;

	//now parallel edges are consecutive in sorted_edges
	bool firstEdge = true;
	for (const auto &ei : sorted_edges)
	{//for
		edge e_act = ei.get_edge();
		int act_s_index = e_act->source()->index();
		int act_t_index = e_act->target()->index();

		if(!firstEdge)
		{//if
			if( (act_s_index == save_s_index && act_t_index == save_t_index) ||
				(act_s_index == save_t_index && act_t_index == save_s_index) )
			{
				if(counter == 1) //first parallel edge
				{
					S.pushBack(e_save);
					new_edgelength[e_save] = E[original_edge[e_save]].get_length() +
						E[original_edge[e_act]].get_length();
				}
				else //more then two parallel edges
					new_edgelength[e_save] +=E[original_edge[e_act]].get_length();

				E[original_edge[e_act]].set_copy_edge(nullptr); //mark copy of edge as deleted
				G_reduced.delEdge(e_act);                    //delete copy edge in G_reduced
				counter++;
			}
			else
			{
				if (counter > 1)
				{
					new_edgelength[e_save]/=counter;
					counter = 1;
				}
				save_s_index = act_s_index;
				save_t_index = act_t_index;
				e_save = e_act;
			}
		}//if
		else //first edge
		{
			firstEdge = false;
			save_s_index = act_s_index;
			save_t_index = act_t_index;
			e_save = e_act;
		}
	}//for

	//treat special case (last edges were multiple edges)
	if(counter >1)
		new_edgelength[e_save]/=counter;
}


void FMMMLayout::update_edgelength(
	List<edge>& S,
	EdgeArray<double>& new_edgelength,
	EdgeArray<EdgeAttributes>& E_reduced)
{
	while (!S.empty())
	{
		edge e = S.popFrontRet();
		E_reduced[e].set_length(new_edgelength[e]);
	}
}


#if 0
inline double FMMMLayout::get_post_rep_force_strength(int n)
{
	return min(0.2,400.0/double(n));
}
#endif


void FMMMLayout::make_positions_integer(Graph& G, NodeArray<NodeAttributes>& A)
{
	if (allowedPositions() == apInteger)
	{
		//calculate value of max_integer_position
		max_integer_position = 100 * average_ideal_edgelength * G.numberOfNodes() *
			G.numberOfNodes();
	}

	//restrict positions to lie in [-max_integer_position,max_integer_position]
	//X [-max_integer_position,max_integer_position]
	for (node v : G.nodes) {
		if ((A[v].get_x() > max_integer_position) ||
			(A[v].get_y() > max_integer_position) ||
			(A[v].get_x() < max_integer_position * (-1.0)) ||
			(A[v].get_y() < max_integer_position * (-1.0)))
		{
			DPoint cross_point;
			DPoint nullpoint(0, 0);
			DPoint old_pos(A[v].get_x(), A[v].get_y());
			DPoint lt(max_integer_position * (-1.0), max_integer_position);
			DPoint rt(max_integer_position, max_integer_position);
			DPoint lb(max_integer_position * (-1.0), max_integer_position * (-1.0));
			DPoint rb(max_integer_position, max_integer_position * (-1.0));
			DLine s(nullpoint, old_pos);
			DLine left_bound(lb, lt);
			DLine right_bound(rb, rt);
			DLine top_bound(lt, rt);
			DLine bottom_bound(lb, rb);

			if (s.intersection(left_bound, cross_point))
			{
				A[v].set_x(cross_point.m_x);
				A[v].set_y(cross_point.m_y);
			}
			else if (s.intersection(right_bound, cross_point))
			{
				A[v].set_x(cross_point.m_x);
				A[v].set_y(cross_point.m_y);
			}
			else if (s.intersection(top_bound, cross_point))
			{
				A[v].set_x(cross_point.m_x);
				A[v].set_y(cross_point.m_y);
			}
			else if (s.intersection(bottom_bound, cross_point))
			{
				A[v].set_x(cross_point.m_x);
				A[v].set_y(cross_point.m_y);
			}
			else cout << "Error FMMMLayout:: make_positions_integer()" << endl;
		}
	}
	//make positions integer
	for (node v : G.nodes)
	{
		double new_x = floor(A[v].get_x());
		double new_y = floor(A[v].get_y());
		if (new_x < down_left_corner.m_x)
		{
			boxlength += 2;
			down_left_corner.m_x = down_left_corner.m_x - 2;
		}
		if (new_y < down_left_corner.m_y)
		{
			boxlength += 2;
			down_left_corner.m_y = down_left_corner.m_y - 2;
		}
		A[v].set_x(new_x);
		A[v].set_y(new_y);
	}
}


void FMMMLayout::create_postscript_drawing(GraphAttributes& AG, char* ps_file)
{
	ofstream out_fmmm(ps_file, ios::out);
	if (!ps_file)
		cout << ps_file << " could not be opened !" << endl;

	const Graph& G = AG.constGraph();
	double x_min = AG.x(G.firstNode());
	double x_max = x_min;
	double y_min = AG.y(G.firstNode());
	double y_max = y_min;
	double max_dist;
	double scale_factor;

	for (node v : G.nodes)
	{
		if (AG.x(v) < x_min)
			x_min = AG.x(v);
		else if (AG.x(v) > x_max)
			x_max = AG.x(v);
		if (AG.y(v) < y_min)
			y_min = AG.y(v);
		else if (AG.y(v) > y_max)
			y_max = AG.y(v);
	}
	max_dist = max(x_max - x_min, y_max - y_min);
	scale_factor = 500.0 / max_dist;

	out_fmmm << "%!PS-Adobe-2.0 " << endl;
	out_fmmm << "%%Pages:  1 " << endl;
	out_fmmm << "% %BoundingBox: " << x_min << " " << x_max << " " << y_min << " " << y_max << endl;
	out_fmmm << "%%EndComments " << endl;
	out_fmmm << "%%" << endl;
	out_fmmm << "%% Circle" << endl;
	out_fmmm << "/ellipse_dict 4 dict def" << endl;
	out_fmmm << "/ellipse {" << endl;
	out_fmmm << "  ellipse_dict" << endl;
	out_fmmm << "  begin" << endl;
	out_fmmm << "   newpath" << endl;
	out_fmmm << "   /yrad exch def /xrad exch def /ypos exch def /xpos exch def" << endl;
	out_fmmm << "   matrix currentmatrix" << endl;
	out_fmmm << "   xpos ypos translate" << endl;
	out_fmmm << "   xrad yrad scale" << endl;
	out_fmmm << "  0 0 1 0 360 arc" << endl;
	out_fmmm << "  setmatrix" << endl;
	out_fmmm << "  closepath" << endl;
	out_fmmm << " end" << endl;
	out_fmmm << "} def" << endl;
	out_fmmm << "%% Nodes" << endl;
	out_fmmm << "/v { " << endl;
	out_fmmm << " /y exch def" << endl;
	out_fmmm << " /x exch def" << endl;
	out_fmmm << "1.000 1.000 0.894 setrgbcolor" << endl;
	out_fmmm << "x y 10.0 10.0 ellipse fill" << endl;
	out_fmmm << "0.000 0.000 0.000 setrgbcolor" << endl;
	out_fmmm << "x y 10.0 10.0 ellipse stroke" << endl;
	out_fmmm << "} def" << endl;
	out_fmmm << "%% Edges" << endl;
	out_fmmm << "/e { " << endl;
	out_fmmm << " /b exch def" << endl;
	out_fmmm << " /a exch def" << endl;
	out_fmmm << " /y exch def" << endl;
	out_fmmm << " /x exch def" << endl;
	out_fmmm << "x y moveto a b lineto stroke" << endl;
	out_fmmm << "} def" << endl;
	out_fmmm << "%% " << endl;
	out_fmmm << "%% INIT " << endl;
	out_fmmm << "20  200 translate" << endl;
	out_fmmm << scale_factor << "  " << scale_factor << "  scale " << endl;
	out_fmmm << "1 setlinewidth " << endl;
	out_fmmm << "%%BeginProgram " << endl;

	for (edge e : G.edges)
		out_fmmm << AG.x(e->source()) << " " << AG.y(e->source()) << " "
		<< AG.x(e->target()) << " " << AG.y(e->target()) << " e" << endl;

	for (node v : G.nodes)
		out_fmmm << AG.x(v) << " " << AG.y(v) << " v" << endl;

	out_fmmm << "%%EndProgram " << endl;
	out_fmmm << "showpage " << endl;
	out_fmmm << "%%EOF " << endl;
}


//------------------------- functions for divide et impera step -----------------------

void FMMMLayout::create_maximum_connected_subGraphs(
	Graph& G,
	NodeArray<NodeAttributes>& A,
	EdgeArray<EdgeAttributes>&E,
	Graph G_sub [],
	NodeArray<NodeAttributes> A_sub [],
	EdgeArray<EdgeAttributes> E_sub [],
	NodeArray<int>& component)
{
	//create the subgraphs and save links to subgraph nodes/edges in A
	for (node v_orig : G.nodes)
		A[v_orig].set_subgraph_node(G_sub[component[v_orig]].newNode());

	for (edge e_orig : G.edges)
	{
		node u_orig = e_orig->source();
		node v_orig = e_orig->target();
		E[e_orig].set_subgraph_edge(G_sub[component[u_orig]].newEdge
			(A[u_orig].get_subgraph_node(), A[v_orig].get_subgraph_node()));
	}

	//make A_sub,E_sub valid for the subgraphs
	for (int i = 0; i < number_of_components; i++)
	{
		A_sub[i].init(G_sub[i]);
		E_sub[i].init(G_sub[i]);
	}

	//import information for A_sub,E_sub and links to the original nodes/edges
	//of the subGraph nodes/edges

	for (node v_orig : G.nodes)
	{
		node v_sub = A[v_orig].get_subgraph_node();
		A_sub[component[v_orig]][v_sub].set_NodeAttributes(A[v_orig].get_width(),
			A[v_orig].get_height(), A[v_orig].get_position(),
			v_orig, nullptr);
	}
	for (edge e_orig : G.edges)
	{
		edge e_sub = E[e_orig].get_subgraph_edge();
		node v_orig = e_orig->source();
		E_sub[component[v_orig]][e_sub].set_EdgeAttributes(E[e_orig].get_length(), e_orig, nullptr);
	}
}


void FMMMLayout::pack_subGraph_drawings(
	NodeArray<NodeAttributes>& A,
	Graph G_sub[],
	NodeArray<NodeAttributes> A_sub[])
{
	double aspect_ratio_area, bounding_rectangles_area;
	energybased::MAARPacking P;
	List<Rectangle> R;

	if(stepsForRotatingComponents() == 0) //no rotation
		calculate_bounding_rectangles_of_components(R,G_sub,A_sub);
	else
		rotate_components_and_calculate_bounding_rectangles(R,G_sub,A_sub);

	P.pack_rectangles_using_Best_Fit_strategy(R,pageRatio(),presortCCs(),
		tipOverCCs(),aspect_ratio_area,
		bounding_rectangles_area);
	export_node_positions(A,R,G_sub,A_sub);
}


void FMMMLayout::calculate_bounding_rectangles_of_components(
	List<Rectangle>& R,
	Graph G_sub[],
	NodeArray<NodeAttributes> A_sub[])
{
	int i;
	Rectangle r;
	R.clear();

	for(i=0;i<number_of_components;i++)
	{
		r = calculate_bounding_rectangle(G_sub[i],A_sub[i],i);
		R.pushBack(r);
	}
}


energybased::Rectangle FMMMLayout::calculate_bounding_rectangle(
	Graph& G,
	NodeArray<NodeAttributes>& A,
	int componenet_index)
{
	Rectangle r;
	node v = G.firstNode();
	// max_boundary is the maximum of half of the width and half of the
	// height of each node; (needed to be able to tip rectangles over
	// without having access to the height and width of each node)
	double max_boundary = max(A[v].get_width() / 2, A[v].get_height() / 2);
	double
	  x_min = A[v].get_x() - max_boundary,
	  x_max = A[v].get_x() + max_boundary,
	  y_min = A[v].get_y() - max_boundary,
	  y_max = A[v].get_y() + max_boundary;

	for (v = v->succ(); v; v = v->succ()) {
		max_boundary = max(A[v].get_width()/2, A[v].get_height()/2);
		const double
		  act_x_min = A[v].get_x() - max_boundary,
		  act_x_max = A[v].get_x() + max_boundary,
		  act_y_min = A[v].get_y() - max_boundary,
		  act_y_max = A[v].get_y() + max_boundary;
		if (act_x_min < x_min) x_min = act_x_min;
		if (act_x_max > x_max) x_max = act_x_max;
		if (act_y_min < y_min) y_min = act_y_min;
		if (act_y_max > y_max) y_max = act_y_max;
	}

	//add offset
	x_min -= minDistCC()/2;
	x_max += minDistCC()/2;
	y_min -= minDistCC()/2;
	y_max += minDistCC()/2;

	r.set_rectangle(x_max-x_min,y_max-y_min,x_min,y_min,componenet_index);
	return r;
}


void FMMMLayout::rotate_components_and_calculate_bounding_rectangles(
	List<Rectangle>&R,
	Graph G_sub [],
	NodeArray<NodeAttributes> A_sub [])
{
	Array<NodeArray<DPoint> > best_coords(number_of_components);
	Array<NodeArray<DPoint> > old_coords(number_of_components);
#if 0
	node v_sub;
#endif
	Rectangle r_act, r_best;
	DPoint new_pos, new_dlc;

	R.clear(); //make R empty

	for (int i = 0; i < number_of_components; i++)
	{//allcomponents

		//init r_best, best_area and best_(old)coords
		r_best = calculate_bounding_rectangle(G_sub[i], A_sub[i], i);
		double best_area = calculate_area(r_best.get_width(), r_best.get_height(),
			number_of_components);
		best_coords[i].init(G_sub[i]);
		old_coords[i].init(G_sub[i]);

		for (node v_sub : G_sub[i].nodes)
			old_coords[i][v_sub] = best_coords[i][v_sub] = A_sub[i][v_sub].get_position();

		//rotate the components
		for (int j = 1; j <= stepsForRotatingComponents(); j++)
		{
			//calculate new positions for the nodes, the new rectangle and area
			double angle = Math::pi_2 * (double(j) / double(stepsForRotatingComponents() + 1));
			double sin_j = sin(angle);
			double cos_j = cos(angle);
			for (node v_sub : G_sub[i].nodes)
			{
				new_pos.m_x = cos_j * old_coords[i][v_sub].m_x
					- sin_j * old_coords[i][v_sub].m_y;
				new_pos.m_y = sin_j * old_coords[i][v_sub].m_x
					+ cos_j * old_coords[i][v_sub].m_y;
				A_sub[i][v_sub].set_position(new_pos);
			}

			r_act = calculate_bounding_rectangle(G_sub[i], A_sub[i], i);
			double act_area = calculate_area(r_act.get_width(), r_act.get_height(),
				number_of_components);

			double act_area_PI_half_rotated;
			if (number_of_components == 1)
				act_area_PI_half_rotated = calculate_area(r_act.get_height(),
				r_act.get_width(),
				number_of_components);

			//store placement of the nodes with minimal area (in case that
			//number_of_components >1) else store placement with minimal aspect ratio area
			if (act_area < best_area)
			{
				r_best = r_act;
				best_area = act_area;
				for (node v_sub : G_sub[i].nodes)
					best_coords[i][v_sub] = A_sub[i][v_sub].get_position();
			}
			else if ((number_of_components == 1) && (act_area_PI_half_rotated < best_area))
			{ //test if rotating further with PI_half would be an improvement
				r_best = r_act;
				best_area = act_area_PI_half_rotated;
				for (node v_sub : G_sub[i].nodes)
					best_coords[i][v_sub] = A_sub[i][v_sub].get_position();
				//the needed rotation step follows in the next if statement
			}
		}

		//tipp the smallest rectangle over by angle PI/2 around the origin if it makes the
		//aspect_ratio of r_best more similar to the desired aspect_ratio
		double ratio = r_best.get_width() / r_best.get_height();

		if ((pageRatio() <  1 && ratio > 1) || (pageRatio() >= 1 && ratio < 1))
		{
			for (node v_sub : G_sub[i].nodes)
			{
				new_pos.m_x = best_coords[i][v_sub].m_y*(-1);
				new_pos.m_y = best_coords[i][v_sub].m_x;
				best_coords[i][v_sub] = new_pos;
			}

			//calculate new rectangle
			new_dlc.m_x = r_best.get_old_dlc_position().m_y*(-1) - r_best.get_height();
			new_dlc.m_y = r_best.get_old_dlc_position().m_x;

			double new_width = r_best.get_height();
			double new_height = r_best.get_width();
			r_best.set_width(new_width);
			r_best.set_height(new_height);
			r_best.set_old_dlc_position(new_dlc);
		}

		//save the computed information in A_sub and R
		for (node v_sub : G_sub[i].nodes)
			A_sub[i][v_sub].set_position(best_coords[i][v_sub]);
		R.pushBack(r_best);

	}//allcomponents
}


void FMMMLayout::export_node_positions(
	NodeArray<NodeAttributes>& A,
	List<Rectangle>&  R,
	Graph G_sub [],
	NodeArray<NodeAttributes> A_sub [])
{
	for (const Rectangle &r : R)
	{
		int i = r.get_component_index();
		if (r.is_tipped_over())
		{
			//calculate tipped coordinates of the nodes
			for (node v_sub : G_sub[i].nodes)
			{
				DPoint tipped_pos(-A_sub[i][v_sub].get_y(), A_sub[i][v_sub].get_x());
				A_sub[i][v_sub].set_position(tipped_pos);
			}
		}

		for (node v_sub : G_sub[i].nodes)
		{
			DPoint newpos = A_sub[i][v_sub].get_position() + r.get_new_dlc_position()
				- r.get_old_dlc_position();
			A[A_sub[i][v_sub].get_original_node()].set_position(newpos);
		}
	}
}


//----------------------- functions for multilevel step -----------------------------

inline int FMMMLayout::get_max_mult_iter(int act_level, int max_level, int node_nr)
{
	int iter;
	if(maxIterChange() == micConstant) //nothing to do
		iter =  fixedIterations();
	else if (maxIterChange() == micLinearlyDecreasing) //linearly decreasing values
	{
		if(max_level == 0)
			iter = fixedIterations() +  ((maxIterFactor()-1) * fixedIterations());
		else
			iter = fixedIterations() + int((double(act_level)/double(max_level) ) *
			(maxIterFactor() - 1) * fixedIterations());
	}
	else //maxIterChange == micRapidlyDecreasing (rapidly decreasing values)
	{
		if(act_level == max_level)
			iter = fixedIterations() + int( (maxIterFactor()-1) * fixedIterations());
		else if(act_level == max_level - 1)
			iter = fixedIterations() + int(0.5 * (maxIterFactor()-1) * fixedIterations());
		else if(act_level == max_level - 2)
			iter = fixedIterations() + int(0.25 * (maxIterFactor()-1) * fixedIterations());
		else //act_level >= max_level - 3
			iter = fixedIterations();
	}

	//helps to get good drawings for small graphs and graphs with few multilevels
	if((node_nr <= 500) && (iter < 100))
		return 100;
	else
		return iter;
}


//-------------------------- functions for force calculation ---------------------------

inline void FMMMLayout::calculate_forces(
	Graph& G,
	NodeArray<NodeAttributes>& A,
	EdgeArray<EdgeAttributes>& E,
	NodeArray<DPoint>& F,
	NodeArray<DPoint>& F_attr,
	NodeArray<DPoint>& F_rep,
	NodeArray<DPoint>& last_node_movement,
	int iter,
	int fine_tuning_step)
{
	if(allowedPositions() != apAll)
		make_positions_integer(G,A);
	calculate_attractive_forces(G,A,E,F_attr);
	calculate_repulsive_forces(G,A,F_rep);
	add_attr_rep_forces(G,F_attr,F_rep,F,iter,fine_tuning_step);
	prevent_oscilations(G,F,last_node_movement,iter);
	move_nodes(G,A,F);
	update_boxlength_and_cornercoordinate(G,A);
}


void FMMMLayout::init_boxlength_and_cornercoordinate(
	Graph& G,
	NodeArray<NodeAttributes>& A)
{
	//boxlength is set

	const double MIN_NODE_SIZE = 10;
	const double BOX_SCALING_FACTOR = 1.1;

	double w = 0, h = 0;
	for (node v : G.nodes)
	{
		w += max(A[v].get_width(), MIN_NODE_SIZE);
		h += max(A[v].get_height(), MIN_NODE_SIZE);
	}

	boxlength = ceil(max(w, h) * BOX_SCALING_FACTOR);

	//down left corner of comp. box is the origin
	down_left_corner.m_x = 0;
	down_left_corner.m_y = 0;
}


void FMMMLayout::create_initial_placement(Graph& G, NodeArray<NodeAttributes>& A)
{
	const int BILLION = 1000000000;

	if (initialPlacementForces() == ipfKeepPositions) // don't change anything
	{
		init_boxlength_and_cornercoordinate(G, A);

	}
	else if (initialPlacementForces() == ipfUniformGrid) //set nodes to the midpoints of a  grid
	{ //(uniform on a grid)
		init_boxlength_and_cornercoordinate(G, A);
		int level = static_cast<int>(ceil(Math::log4(G.numberOfNodes())));
		OGDF_ASSERT(level < 31);
		int m = (1 << level) - 1;
		bool finished = false;
		double blall = boxlength / (m + 1); //boxlength for boxes at the lowest level (depth)
		Array<node> all_nodes(G.numberOfNodes());

		int k = 0;
		for (node v : G.nodes)
		{
			all_nodes[k] = v;
			k++;
		}
		node v = all_nodes[0];
		k = 0;
		int i = 0;
		while ((!finished) && (i <= m))
		{
			int j = 0;
			while ((!finished) && (j <= m))
			{
				A[v].set_x(boxlength*i / (m + 1) + blall / 2);
				A[v].set_y(boxlength*j / (m + 1) + blall / 2);
				if (k == G.numberOfNodes() - 1)
					finished = true;
				else
				{
					k++;
					v = all_nodes[k];
				}
				j++;
			}
			i++;
		}
	} //(uniform on a grid)
	else //randomised distribution of the nodes;
	{ //(random)
		init_boxlength_and_cornercoordinate(G, A);
		if (initialPlacementForces() == ipfRandomTime)//(RANDOM based on actual CPU-time)
			srand((unsigned int) time(nullptr));
		else if (initialPlacementForces() == ipfRandomRandIterNr)//(RANDOM based on seed)
			srand(randSeed());

		for (node v : G.nodes)
		{
			DPoint rndp;
			rndp.m_x = double(randomNumber(0, BILLION)) / BILLION; //rand_x in [0,1]
			rndp.m_y = double(randomNumber(0, BILLION)) / BILLION; //rand_y in [0,1]
			A[v].set_x(rndp.m_x*(boxlength - 2) + 1);
			A[v].set_y(rndp.m_y*(boxlength - 2) + 1);
		}
	} //(random)
	update_boxlength_and_cornercoordinate(G, A);
}


void FMMMLayout::init_F(Graph& G, NodeArray<DPoint>& F)
{
	DPoint nullpoint(0, 0);
	for (node v : G.nodes)
		F[v] = nullpoint;
}


void FMMMLayout::make_initialisations_for_rep_calc_classes(Graph& G)
{
	if (repulsiveForcesCalculation() == rfcExact)
		FR.make_initialisations(boxlength, down_left_corner, frGridQuotient());
	else if (repulsiveForcesCalculation() == rfcGridApproximation)
		FR.make_initialisations(boxlength, down_left_corner, frGridQuotient());
	else //(repulsiveForcesCalculation() == rfcNMM
		NM.make_initialisations(G, boxlength, down_left_corner,
		nmParticlesInLeaves(), nmPrecision(),
		nmTreeConstruction(), nmSmallCell());
}


void FMMMLayout::calculate_attractive_forces(
	Graph& G,
	NodeArray<NodeAttributes> & A,
	EdgeArray<EdgeAttributes> & E,
	NodeArray<DPoint>& F_attr)
{
	energybased::numexcept N;
	DPoint f_u;
	DPoint nullpoint (0,0);

	//initialisation
	init_F(G,F_attr);

	//calculation
	for(edge e : G.edges)
	{
		node u = e->source();
		node v = e->target();
		DPoint vector_v_minus_u  = A[v].get_position() - A[u].get_position();
		double norm_v_minus_u = vector_v_minus_u.norm();
		if(vector_v_minus_u == nullpoint)
			f_u = nullpoint;
		else if(!N.f_near_machine_precision(norm_v_minus_u,f_u))
		{
			double scalar = f_attr_scalar(norm_v_minus_u,E[e].get_length())/norm_v_minus_u;
			f_u.m_x = scalar * vector_v_minus_u.m_x;
			f_u.m_y = scalar * vector_v_minus_u.m_y;
		}

		F_attr[v] = F_attr[v] - f_u;
		F_attr[u] = F_attr[u] + f_u;
	}
}


double FMMMLayout::f_attr_scalar(double d, double ind_ideal_edge_length)
{
	double s(0);

	switch (forceModel()) {
	case fmFruchtermanReingold:
		s =  d*d/(ind_ideal_edge_length*ind_ideal_edge_length*ind_ideal_edge_length);
		break;
	case fmEades:
		{
			const double c = 10;
			if (d == 0)
				s = -1e10;
			else
				s =  c * std::log2(d/ind_ideal_edge_length) / ind_ideal_edge_length;
			break;
		}
	case fmNew:
		{
			const double c =  std::log2(d/ind_ideal_edge_length);
			if (d > 0)
				s =  c * d * d /
				(ind_ideal_edge_length * ind_ideal_edge_length * ind_ideal_edge_length);
			else
				s = -1e10;
			break;
		}
	default:
		cerr << "Error FMMMLayout::f_attr_scalar" << endl;
		OGDF_ASSERT(false);
	}

	return s;
}


void FMMMLayout::add_attr_rep_forces(
	Graph& G,
	NodeArray<DPoint>& F_attr,
	NodeArray<DPoint>& F_rep,
	NodeArray<DPoint>& F,
	int iter,
	int fine_tuning_step)
{
	energybased::numexcept N;
	DPoint nullpoint(0, 0);

	//set cool_factor
	if (coolTemperature() == false)
		cool_factor = 1.0;
	else if ((coolTemperature() == true) && (fine_tuning_step == 0))
	{
		if (iter == 1)
			cool_factor = coolValue();
		else
			cool_factor *= coolValue();
	}

	if (fine_tuning_step == 1)
		cool_factor /= 10.0; //decrease the temperature rapidly
	else if (fine_tuning_step == 2)
	{
		if (iter <= fineTuningIterations() - 5)
			cool_factor = fineTuneScalar(); //decrease the temperature rapidly
		else
			cool_factor = (fineTuneScalar() / 10.0);
	}

	//set the values for the spring strength and strength of the rep. force field
	double act_spring_strength, act_rep_force_strength;
	if (fine_tuning_step <= 1)//usual case
	{
		act_spring_strength = springStrength();
		act_rep_force_strength = repForcesStrength();
	}
	else if (!adjustPostRepStrengthDynamically())
	{
		act_spring_strength = postSpringStrength();
		act_rep_force_strength = postStrengthOfRepForces();
	}
	else //adjustPostRepStrengthDynamically())
	{
		act_spring_strength = postSpringStrength();
		act_rep_force_strength = get_post_rep_force_strength(G.numberOfNodes());
	}

	for (node v : G.nodes)
	{
		DPoint f;
		f.m_x = act_spring_strength * F_attr[v].m_x + act_rep_force_strength * F_rep[v].m_x;
		f.m_y = act_spring_strength * F_attr[v].m_y + act_rep_force_strength * F_rep[v].m_y;
		f.m_x = average_ideal_edgelength * average_ideal_edgelength * f.m_x;
		f.m_y = average_ideal_edgelength * average_ideal_edgelength * f.m_y;

		double norm_f = f.norm();

		DPoint force;
		if (f == nullpoint)
			force = nullpoint;
		else if (N.f_near_machine_precision(norm_f, force))
			restrict_force_to_comp_box(force);
		else
		{
			double scalar = min(norm_f * cool_factor * forceScalingFactor(), max_radius(iter)) / norm_f;
			force.m_x = scalar * f.m_x;
			force.m_y = scalar * f.m_y;
		}
		F[v] = force;
	}
}


void FMMMLayout::move_nodes(
	Graph& G,
	NodeArray<NodeAttributes>& A,
	NodeArray<DPoint>& F)
{
	for(node v : G.nodes)
		A[v].set_position(A[v].get_position() + F[v]);
}


void FMMMLayout::update_boxlength_and_cornercoordinate(
	Graph& G,
	NodeArray<NodeAttributes>&A)
{
	node vFirst = G.firstNode();
	DPoint midpoint = A[vFirst].get_position();

	double xmin, xmax, ymin, ymax;
	xmin = xmax = midpoint.m_x;
	ymin = ymax = midpoint.m_y;

	for (node v : G.nodes)
	{
		midpoint = A[v].get_position();
		if (midpoint.m_x < xmin)
			xmin = midpoint.m_x;
		if (midpoint.m_x > xmax)
			xmax = midpoint.m_x;
		if (midpoint.m_y < ymin)
			ymin = midpoint.m_y;
		if (midpoint.m_y > ymax)
			ymax = midpoint.m_y;
	}

	//set down_left_corner and boxlength

	down_left_corner.m_x = floor(xmin - 1);
	down_left_corner.m_y = floor(ymin - 1);
	boxlength = ceil(max(ymax - ymin, xmax - xmin) *1.01 + 2);

	//exception handling: all nodes have same x and y coordinate
	if (boxlength <= 2)
	{
		boxlength = G.numberOfNodes() * 20;
		down_left_corner.m_x = floor(xmin) - (boxlength / 2);
		down_left_corner.m_y = floor(ymin) - (boxlength / 2);
	}

	//export the boxlength and down_left_corner values to the rep. calc. classes

	if (repulsiveForcesCalculation() == rfcExact ||
		repulsiveForcesCalculation() == rfcGridApproximation)
		FR.update_boxlength_and_cornercoordinate(boxlength, down_left_corner);
	else //repulsiveForcesCalculation() == rfcNMM
		NM.update_boxlength_and_cornercoordinate(boxlength, down_left_corner);
}


void FMMMLayout::set_average_ideal_edgelength(
	Graph& G,
	EdgeArray<EdgeAttributes>& E)
{
	if(G.numberOfEdges() > 0)
	{
		double averagelength = 0;
		for(edge e : G.edges)
			averagelength += E[e].get_length();
		average_ideal_edgelength = averagelength/G.numberOfEdges();
	}
	else
		average_ideal_edgelength = 50;
}


double FMMMLayout::get_average_forcevector_length(Graph& G, NodeArray<DPoint>& F)
{
	double lengthsum = 0;
	for (node v : G.nodes)
		lengthsum += F[v].norm();
	lengthsum /= G.numberOfNodes();
	return lengthsum;
}


void FMMMLayout::prevent_oscilations(
	Graph& G,
	NodeArray<DPoint>& F,
	NodeArray<DPoint>& last_node_movement,
	int iter)
{
	const double pi_times_1_over_6 = 0.52359878;
	const double pi_times_2_over_6 = 2 * pi_times_1_over_6;
	const double pi_times_3_over_6 = 3 * pi_times_1_over_6;
	const double pi_times_4_over_6 = 4 * pi_times_1_over_6;
	const double pi_times_5_over_6 = 5 * pi_times_1_over_6;
	const double pi_times_7_over_6 = 7 * pi_times_1_over_6;
	const double pi_times_8_over_6 = 8 * pi_times_1_over_6;
	const double pi_times_9_over_6 = 9 * pi_times_1_over_6;
	const double pi_times_10_over_6 = 10 * pi_times_1_over_6;
	const double pi_times_11_over_6 = 11 * pi_times_1_over_6;

	DPoint nullpoint (0,0);

	if (iter > 1) // usual case
	{
		for(node v : G.nodes)
		{
			DPoint force_new (F[v].m_x,F[v].m_y);
			DPoint force_old (last_node_movement[v].m_x,last_node_movement[v].m_y);
			double norm_new = F[v].norm();
			double norm_old  = last_node_movement[v].norm();
			if ((norm_new > 0) && (norm_old > 0))
			{
				double quot_old_new =  norm_old / norm_new;

				// prevent oszilations
				// angle in [0,2pi) measured counterclockwise
				double fi = angle(nullpoint,force_old,force_new);
				if( ( fi <= pi_times_1_over_6 || fi >= pi_times_11_over_6 ) && norm_new > norm_old*2.0 )
				{
					F[v].m_x = quot_old_new * 2.0 * F[v].m_x;
					F[v].m_y = quot_old_new * 2.0 * F[v].m_y;
				}
				else if ( fi >= pi_times_1_over_6 && fi <= pi_times_2_over_6 && norm_new > norm_old*1.5 )
				{
					F[v].m_x = quot_old_new * 1.5 * F[v].m_x;
					F[v].m_y = quot_old_new * 1.5 * F[v].m_y;
				}
				else if ( fi >= pi_times_2_over_6 && fi <= pi_times_3_over_6 && norm_new > norm_old )
				{
					F[v].m_x = quot_old_new * F[v].m_x;
					F[v].m_y = quot_old_new * F[v].m_y;
				}
				else if ( fi >= pi_times_3_over_6 && fi <= pi_times_4_over_6 && norm_new > norm_old*0.66666666 )
				{
					F[v].m_x = quot_old_new * 0.66666666 * F[v].m_x;
					F[v].m_y = quot_old_new * 0.66666666 * F[v].m_y;
				}
				else if ( fi >= pi_times_4_over_6 && fi <= pi_times_5_over_6 && norm_new > norm_old*0.5 )
				{
					F[v].m_x = quot_old_new * 0.5 * F[v].m_x;
					F[v].m_y = quot_old_new * 0.5 * F[v].m_y;
				}
				else if ( fi >= pi_times_5_over_6 && fi <= pi_times_7_over_6 && norm_new > norm_old*0.33333333 )
				{
					F[v].m_x = quot_old_new * 0.33333333 * F[v].m_x;
					F[v].m_y = quot_old_new * 0.33333333 * F[v].m_y;
				}
				else if ( fi >= pi_times_7_over_6 && fi <= pi_times_8_over_6 && norm_new > norm_old*0.5 )
				{
					F[v].m_x = quot_old_new * 0.5 * F[v].m_x;
					F[v].m_y = quot_old_new * 0.5 * F[v].m_y;
				}
				else if ( fi >= pi_times_8_over_6 && fi <= pi_times_9_over_6 && norm_new > norm_old*0.66666666 )
				{
					F[v].m_x = quot_old_new * 0.66666666 * F[v].m_x;
					F[v].m_y = quot_old_new * 0.66666666 * F[v].m_y;
				}
				else if ( fi >= pi_times_9_over_6 && fi <= pi_times_10_over_6 && norm_new > norm_old )
				{
					F[v].m_x = quot_old_new * F[v].m_x;
					F[v].m_y = quot_old_new * F[v].m_y;
				}
				else if ( fi >= pi_times_10_over_6 && fi <= pi_times_11_over_6 && norm_new > norm_old*1.5 )
				{
					F[v].m_x = quot_old_new * 1.5 * F[v].m_x;
					F[v].m_y = quot_old_new * 1.5 * F[v].m_y;
				}
			}
			last_node_movement[v]= F[v];
		}
	}
	else if (iter == 1)
		init_last_node_movement(G,F,last_node_movement);
}


double FMMMLayout::angle(DPoint& P, DPoint& Q, DPoint& R)
{
	double dx1 = Q.m_x - P.m_x;
	double dy1 = Q.m_y - P.m_y;
	double dx2 = R.m_x - P.m_x;
	double dy2 = R.m_y - P.m_y;
	double fi;//the angle

	if ((dx1 == 0 && dy1 == 0) || (dx2 == 0 && dy2 == 0))
		cout<<"Multilevel::angle()"<<endl;

	double norm  = (dx1*dx1+dy1*dy1)*(dx2*dx2+dy2*dy2);
	double cosfi = (dx1*dx2+dy1*dy2) / sqrt(norm);

	if (cosfi <= -1.0 ) fi = Math::pi;
	else
	{
		fi = acos(cosfi);
		if (dx1*dy2 < dy1*dx2) fi = -fi;
		if (fi < 0) fi += 2*Math::pi;
	}
	return fi;
}


void FMMMLayout::init_last_node_movement(
	Graph& G,
	NodeArray<DPoint>& F,
	NodeArray<DPoint>& last_node_movement)
{
	for(node v : G.nodes)
		last_node_movement[v]= F[v];
}


void FMMMLayout::adapt_drawing_to_ideal_average_edgelength(
	Graph& G,
	NodeArray<NodeAttributes>& A,
	EdgeArray<EdgeAttributes>& E)
{
	double sum_real_edgelength = 0;
	double sum_ideal_edgelength = 0;
	for (edge e : G.edges)
	{
		sum_ideal_edgelength += E[e].get_length();
		sum_real_edgelength += (A[e->source()].get_position() - A[e->target()].get_position()).norm();
	}

	double area_scaling_factor;
	if (sum_real_edgelength == 0) // very very unlikly case
		area_scaling_factor = 1;
	else
		area_scaling_factor = sum_ideal_edgelength / sum_real_edgelength;

	DPoint new_pos;
	for (node v : G.nodes)
	{
		new_pos.m_x = resizingScalar() * area_scaling_factor * A[v].get_position().m_x;
		new_pos.m_y = resizingScalar() * area_scaling_factor * A[v].get_position().m_y;
		A[v].set_position(new_pos);
	}
}


} //end namespace ogdf
