/*!
 * 
 *
 * \brief       Calculates the hypervolume covered by a front of non-dominated points.
 * 
 * 
 *
 * \author      T.Voss
 * \date        2010
 *
 *
 * \par Copyright 1995-2017 Shark Development Team
 * 
 * <BR><HR>
 * This file is part of Shark.
 * <http://shark-ml.org/>
 * 
 * Shark is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Shark is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with Shark.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef SHARK_ALGORITHMS_DIRECTSEARCH_INDICATORS_HYPERVOLUMEINDICATOR_H
#define SHARK_ALGORITHMS_DIRECTSEARCH_INDICATORS_HYPERVOLUMEINDICATOR_H

#include "HypervolumeContribution.h"

#include <algorithm>
#include <vector>
#include <numeric>
using namespace std;
namespace shark {

///  \brief Calculates the hypervolume covered by a front of non-dominated points.
///
/// If given, the Indicator uses a provided reference value that can be set via setReference. 
/// Otherwise, it is computed from the data by using the maximum value in the set. As this usually
/// gives 0 contribution to the extremum points (i.e. the ones with best function value), those
/// points are skipped when computing the contribution (i.e. extremum points are never selected).
/// Note, that for boundary points that are not extrema, this does not hold and they are selected.
///
/// for problems with many objectives, an approximative algorithm can be used.
struct HypervolumeIndicator {
	/// \brief Determines the point contributing the least hypervolume to the overall front of points.
	///
	size_t leastContributor( std::vector<vector<double> > const& front, std::vector<double> const& /*archive*/)const{
		HypervolumeContribution algorithm;
		if(m_reference.size() != 0)
			return m_algorithm.smallest(front,1,m_reference)[0].second;
		else	
			return m_algorithm.smallest(front,1)[0].second;
	}
	
	std::vector<std::size_t> leastContributors( std::vector<vector<double> > const& front, std::vector<vector<double>> const& archive, std::size_t K)const{
		std::vector<std::size_t> indices;
		std::vector<vector<double>> points(front.begin(),front.end());
		std::vector<std::size_t> activeIndices(points.size());
		std::iota(activeIndices.begin(),activeIndices.end(),0);
		for(std::size_t k=0; k != K; ++k){
			std::size_t index = -1;//leastContributor(points,archive);
			points.erase(points.begin()+index);
			indices.push_back(activeIndices[index]);
			activeIndices.erase(activeIndices.begin()+index);
		}
		return indices;
	}
	std::vector<std::pair<double,std::size_t> > leastContributors( std::vector<vector<double> > const& front, std::size_t K)const{
	//std::vector<KeyValuePair<double,std::size_t> > leastContributors( std::vector<double> const& front, std::size_t K)const{
		if(m_reference.size() != 0)
			return m_algorithm.smallest(front,K,m_reference);
		else	
			return m_algorithm.smallest(front,K);
	}
	std::vector<std::pair<double,std::size_t> > lastIdxContributor( std::vector<vector<double> > const& front, std::size_t K)const{
	//std::vector<KeyValuePair<double,std::size_t> > leastContributors( std::vector<double> const& front, std::size_t K)const{
		if(m_reference.size() != 0)
			return m_algorithm.lastIdxContributor(front,K,m_reference);
		else	
			return m_algorithm.lastIdxContributor(front,K);
	}
	void init(std::size_t /*numOfObjectives*/, std::size_t /*mu*/, unsigned int& /*rng*/){}
	
	/// \brief Sets the reference point. 
	///
	/// If no point is set, it is estimated from the current front and the extremum points are never selected.
	void setReference(std::vector<double> const& newReference){
		m_reference = newReference;
	}
	
	/// \brief Whether the approximtive algorithm should be used on large problems
//	void useApproximation(bool useApproximation){
//		//m_algorithm.useApproximation(useApproximation);
//	}
	
	///\brief Error bound for the approximative algorithm
//	double approximationEpsilon()const{
//		return m_algorithm.approximationEpsilon();
//	}
	///\brief Error bound for the approximative algorithm
//	double& approximationEpsilon(){
//		return m_algorithm.approximationEpsilon();
//	}
	
	///\brief Error probability for the approximative algorithm
//	double approximationDelta()const{
//		return m_algorithm.approximationDelta();
//	}
	
	///\brief Error probability for the approximative algorithm
//	double& approximationDelta(){
//		return m_algorithm.approximationDelta();
//	}
//private:
public:
	std::vector<double> m_reference;
	HypervolumeContribution m_algorithm;
};
}

#endif
