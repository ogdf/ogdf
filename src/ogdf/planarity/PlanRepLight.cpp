
#include <ogdf/planarity/PlanRepLight.h>


namespace ogdf {


PlanRepLight::PlanRepLight(const PlanRep &pr)
	: m_ccInfo(pr.ccInfo()), m_pr(pr), m_currentCC(-1), m_eAuxCopy(pr.original())
{
	GraphCopy::createEmpty(pr.original());
}


void PlanRepLight::initCC(int cc)
{
	if (m_currentCC >= 0)
	{
		for(int i = m_ccInfo.startNode(m_currentCC); i < m_ccInfo.stopNode(m_currentCC); ++i)
			m_vCopy[m_ccInfo.v(i)] = 0;

		for(int i = m_ccInfo.startEdge(m_currentCC); i < m_ccInfo.stopEdge(m_currentCC); ++i)
			m_eCopy[m_ccInfo.e(i)].clear();
	}

	m_currentCC = cc;
	initByCC(m_ccInfo, cc, m_eAuxCopy);
}


} // end namespace ogdf
