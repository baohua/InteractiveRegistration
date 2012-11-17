#ifndef __itkGSPolyAffineTransform_hxx
#define __itkGSPolyAffineTransform_hxx

#include "itkGSPolyAffineTransform.h"
#include <iostream>
#include <algorithm>
#include <iterator>
#include <string>
#include "itkBinaryBallStructuringElement.h"
#include "itkSobelEdgeDetectionImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkKdTreeGenerator.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkVector.h"
#include "itkListSample.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkEuclideanDistanceMetric.h"







namespace itk
{
template <class TScalar, unsigned int NDimensions>
GSPolyAffineTransform<TScalar, NDimensions>::GSPolyAffineTransform()
: Superclass()
{

}

template <class TScalar, unsigned int NDimensions>
void GSPolyAffineTransform<TScalar, NDimensions>::ComputeFieldTransfromFromLocalAffineTransform()
{


	SampleFromLabelMask(m_LabelImage, m_SampleSpacing, m_ControlPointList, m_ControlPointLabelList);
	unsigned int nb_pts = m_ControlPointList.size();
//	for(int i=0; i<nb_pts; i++) std::cout << "["<< i << "]" << m_ControlPointList[i] << ":" << m_ControlPointLabelList[i] << std::endl;
	ConvertAffineTransformsToLogarithm();
	// PrecomputeTrajectoryAndCheckCollision();
	PrecomputeTrajectoryAndCheckCollision2();

	// allocate the trajectory mask and helper strucuture for each collsion period;

	unsigned int nbCollision = m_CollisionTimeList.size() - 1;

	m_C.clear();
	m_C.resize(nbCollision);

	for(unsigned int i=0; i<nbCollision; i++) {
	  GetStationaryVieldCopyPasteDecreasing(i, m_C[i].velocityFieldTransform);
	}

}

template <class TScalar, unsigned int NDimensions>
void GSPolyAffineTransform<TScalar, NDimensions>::GetStationaryVieldCopyPasteDecreasing(unsigned int collision_time, ConstantVelocityFieldTransformPointerType &vfield) {

  // most important function here!!
  unsigned int ind_t1 = m_CollisionTimeList[collision_time];
  unsigned int ind_t2 = m_CollisionTimeList[collision_time+1]-1;

  // v = \sum {wi * vi}
  //1. construct the trajectory map between time [ind_t1, ind_t2]
  // typename LabelImageType::Pointer traj ;
  AllocateNewImageOfSameDimension(m_LabelImage, m_C[collision_time].trajectoryImage);
  unsigned int radius = 2; // m_InterTrajectoryDistThres;

  std::cout << "debug!" << std::endl;
  // try using shapedneighborhood iterator
  Size<InputDimension> szradius;
  szradius.Fill(radius);
  SNIteratorType it(szradius, m_C[collision_time].trajectoryImage, m_C[collision_time].trajectoryImage->GetBufferedRegion());
  FillBallShapedNeighborhoodIterator(radius, it);
  std::cout << "debug2!" << std::endl;

  ScatterLabeledPointsIntoTrajectoryImage(collision_time, m_C[collision_time].trajectoryImage, it);

  typename DistanceImageType::Pointer traj_dist_map;
  typename LabelImageType::Pointer traj_dist_map_label;

  PicslImageHelper::WriteResultImage<LabelImageType>(m_C[collision_time].trajectoryImage, "trajectory_mask.nii.gz");
  GetDistanceTransformFromLabelImage(m_C[collision_time].trajectoryImage, m_C[collision_time].trajectoryDistanceImage, m_C[collision_time].trajectoryDistanceVoronoiImage);

  PicslImageHelper::WriteResultImage<DistanceImageType>(m_C[collision_time].trajectoryDistanceImage, "trajectory_mask_distance.nii.gz");
  PicslImageHelper::WriteResultImage<LabelImageType>(m_C[collision_time].trajectoryDistanceVoronoiImage, "trajectory_mask_distance_voronoi.nii.gz");

  unsigned int voronoi_band_radius = ceil(m_InterTrajectoryDistThres / 2);
  GetTrajectoryVoronoiBand(m_C[collision_time].trajectoryDistanceVoronoiImage, voronoi_band_radius, m_C[collision_time].voronoiBandImage, m_C[collision_time].voronoiEdgeImage);

  AssignWeightsToPointsInsideBand(m_C[collision_time].voronoiBandImage, m_C[collision_time].voronoiEdgeImage, collision_time);

  MixingVelocityFieldsFromEveryLabel(collision_time);


  m_C[collision_time].velocityFieldTransform->SetLowerTimeBound(1.0 * ind_t1 / m_NumberOfTimeSteps);
  // +1 here is a must since [0,24] out of 25 steps should map to [0,1]
  m_C[collision_time].velocityFieldTransform->SetUpperTimeBound(1.0 * (ind_t2+1) / m_NumberOfTimeSteps);
  m_C[collision_time].velocityFieldTransform->IntegrateVelocityField();

  std::cout << "NumberOfTimeSteps="<<m_NumberOfTimeSteps << "[" << ind_t1 << "," << ind_t2 << "]" << std::endl;

  typename ConstantVelocityFieldTransformType::DisplacementFieldType::Pointer disp = m_C[collision_time].velocityFieldTransform->GetDisplacementField();
  PicslImageHelper::WriteResultImage<ConstantVelocityFieldType>(disp, "displacement.mhd");


//  typedef itk::Vector<RealType, VImageDimension> VectorType;
//  VectorType zeroVector( 0.0 );
//  typedef itk::Image<VectorType, VImageDimension> DisplacementFieldType;
//  typedef itk::Image<VectorType, VImageDimension> ConstantVelocityFieldType;
//  typename ConstantVelocityFieldType::Pointer displacementField = ConstantVelocityFieldType::New();
//  displacementField->CopyInformation( fixedImage );
//  displacementField->SetRegions( fixedImage->GetBufferedRegion() );
//  displacementField->Allocate();
//  displacementField->FillBuffer( zeroVector );
//
//
//  vfield = ConstantVelocityFieldTransform::New();
//  vfield->SetConstantVelocityField( displacementField );
//
}
template <class TScalar, unsigned int NDimensions>
void GSPolyAffineTransform<TScalar, NDimensions>::MixingVelocityFieldsFromEveryLabel(unsigned int collision_time)
{
  BinaryImagePointerType band = m_C[collision_time].voronoiBandImage;
  LabelImagePointerType voronoi = m_C[collision_time].trajectoryDistanceVoronoiImage;
  LabelImagePointerType labelImage = m_LabelImage;

  m_C[collision_time].velocityFieldTransform = ConstantVelocityFieldTransformType::New();


  typename ConstantVelocityFieldType::Pointer velocityFieldImage = ConstantVelocityFieldType::New();

  ConstantVelocityFieldTransformPointerType velocityFieldTransform = m_C[collision_time].velocityFieldTransform;

  //1. allocate velocityField and velocityFieldTransform
  AllocateNewVectorImageOfSameDimension(labelImage, velocityFieldImage);
  velocityFieldTransform->SetConstantVelocityField(velocityFieldImage);

  //1. iterate over the band image and fill those points outside of band:
  //  v(x) = L_k, k = label
  ImageRegionIteratorWithIndex<ConstantVelocityFieldType> vit(velocityFieldImage, velocityFieldImage->GetBufferedRegion());
  ImageRegionConstIterator<BinaryImageType> bit(band, band->GetBufferedRegion());
  // ImageRegionConstIterator<LabelImageType> lit(labelImage, labelImage->GetBufferedRegion());
  ImageRegionConstIterator<LabelImageType> voit(voronoi, voronoi->GetBufferedRegion());


  unsigned int nbPtsOutOfBand = 0;
  for(vit.GoToBegin(), bit.GoToBegin(), voit.GoToBegin(); !vit.IsAtEnd(); ++vit, ++bit, ++voit) {
    if (bit.Get() == 0) {
      unsigned int label_index = voit.Get() - 1; // label 1 is in the 0-th of m_LocalAffineList
      typename ConstantVelocityFieldType::IndexType index = vit.GetIndex();
      typename ConstantVelocityFieldType::PointType point;
      velocityFieldImage->TransformIndexToPhysicalPoint(index, point);
      VelocityPixelType v = m_LocalAffineList[label_index]->GetVelocityAtPoint(point);
      vit.Set(v);
      ++nbPtsOutOfBand;
    }
  }

  //2. iterate over the bandPointSampleVector, which contains all points inside the band and do averaging here
  // not since there is at least one weight w_k which should be >0.5, so we wont worry about the denominator will be
  // a small number.
  // v = \sum{v_k * w_k} / \sum{w_k}
  unsigned int nbPtsInBand = m_C[collision_time].bandPointSampleVector.size();
  unsigned int nbAffine = m_LocalAffineList.size();
  for(unsigned int i=0; i<nbPtsInBand; i++) {
    typename ConstantVelocityFieldType::IndexType index = m_C[collision_time].bandPointSampleVector[i].index;
    typename ConstantVelocityFieldType::PointType point;
    velocityFieldImage->TransformIndexToPhysicalPoint(index, point);

    unsigned int nbNoneZeroWeights = m_C[collision_time].bandPointSampleVector[i].neighbor_labels.size();
    VelocityPixelType v = NumericTraits<VelocityPixelType>::ZeroValue();

    std::cout << "[" << i << "]";
    for(unsigned int k=0; k<nbNoneZeroWeights; k++) {
      unsigned int label_index = m_C[collision_time].bandPointSampleVector[i].neighbor_labels[k] - 1; //label_index = label - 1, for each affine
      float w = m_C[collision_time].bandPointSampleVector[i].neighbor_weights[k];
      v += m_LocalAffineList[label_index]->GetVelocityAtPoint(point) * w;
      std::cout << m_LocalAffineList[label_index]->GetVelocityAtPoint(point) <<  "*" << w << "+";
    }
    std::cout << " = " << v << std::endl;
    velocityFieldImage->SetPixel(index, v);
  }

  std::cout << "nbPtsOutOfBand=" << nbPtsOutOfBand << " nbPtsInBand="<<nbPtsInBand << " ImageSize" << m_LabelImage->GetBufferedRegion().GetSize() << std::endl;
}

template <class TScalar, unsigned int NDimensions>
void GSPolyAffineTransform<TScalar, NDimensions>::
AssignWeightsToPointsInsideBand(const BinaryImagePointerType &band, const BinaryImagePointerType &edge, unsigned int collision_time)
{
  // for each label
  unsigned int nbAffine = m_LocalAffineList.size();

  // construct kd-tree for points in 'edge' of the same label

    {
    unsigned int nbPtsInEdge = 0;
    ImageRegionConstIterator<BinaryImageType> it( edge, edge->GetBufferedRegion() );
    for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
      {
      if ( it.Get() )
        nbPtsInEdge++;
      }

    std::cout << "points in edge: " << nbPtsInEdge << std::endl;
    }


  typedef typename itk::Vector< float, InputDimension > MeasurementVectorType;
  typedef typename itk::Statistics::ListSample< MeasurementVectorType > SampleType;
  typedef typename itk::Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
  typedef typename TreeGeneratorType::KdTreeType TreeType;
  typedef typename TreeType::NearestNeighbors NeighborsType;
  typedef typename TreeType::KdTreeNodeType NodeType;
  typedef typename TreeType::Pointer TreePointerType;


  typedef itk::Statistics::EuclideanDistanceMetric<MeasurementVectorType>
     DistanceMetricType;
   typedef typename DistanceMetricType::OriginType OriginType;
   typename DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();
   OriginType origin( InputDimension );


  std::vector<TreePointerType> treeListForEachAffine;
  treeListForEachAffine.resize(nbAffine);

  BinaryImagePointerType edgeImage = m_C[collision_time].voronoiEdgeImage;
  LabelImagePointerType voronoiImage = m_C[collision_time].trajectoryDistanceVoronoiImage;
  ImageRegionIteratorWithIndex<BinaryImageType> eit(edgeImage, edgeImage->GetBufferedRegion());
  ImageRegionConstIterator<LabelImageType> lit(voronoiImage, voronoiImage->GetBufferedRegion());

  for(unsigned int i=0; i<nbAffine; i++){
    unsigned int label = i+1;
    // GetKdTreePerLabelAtCollisionTime(collision_time, label, treeListForEachAffine[i]);
    unsigned int nbPtsInEdgeEachAffine = 0;
    for(eit.GoToBegin(), lit.GoToBegin(); !eit.IsAtEnd(); ++eit, ++lit)
      if (eit.Get() && lit.Get()==label) ++nbPtsInEdgeEachAffine;

    typename SampleType::Pointer sample = SampleType::New();
    sample->SetMeasurementVectorSize( InputDimension );
    // sample->Resize(nbPtsInEdgeEachAffine);

    std::cout << "ListSampleSize Before = " << sample->Size() << std::endl;

    MeasurementVectorType mv;
    unsigned int nbPtsInEdgeEachAffine2 = 0;
    for(eit.GoToBegin(), lit.GoToBegin(); !eit.IsAtEnd(); ++eit, ++lit)
      if (eit.Get() && lit.Get()==label) {
        typename LabelImageType::IndexType ind = eit.GetIndex();
        for(unsigned int i=0; i<InputDimension; i++) mv[i] = ind[i];
        sample->PushBack( mv );
        nbPtsInEdgeEachAffine2++;
      }

    std::cout << "ListSampleSize = " << sample->Size() << std::endl;
    std::cout << "nbPtsInEdgeEachAffine2 = " << nbPtsInEdgeEachAffine2 << std::endl;
    std::cout << "nbPtsInEdgeEachAffine = " << nbPtsInEdgeEachAffine << std::endl;

    typename TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
    treeGenerator->SetSample( sample );
    treeGenerator->SetBucketSize( 4 );
    treeGenerator->Update();
    treeListForEachAffine[i] = treeGenerator->GetOutput();

    std::cout << treeListForEachAffine[i] << std::endl;
    typename SampleType::ConstPointer sample1 = treeListForEachAffine[i]->GetSample();
    std::cout << sample1 << std::endl;
    sample1->Register(); //GS: have to be a but in Kdtree, sample will get deleted after moving out.
  }

  // iterate over all points in band and query each kd-tree for closet distance to each label in voronoi

  // first count number of voxels inside a band for allocate bandPointSampleVector
  unsigned int nbPtsInBand = 0;
  ImageRegionConstIterator<BinaryImageType> bit(band, band->GetBufferedRegion());
  for(bit.GoToBegin(); !bit.IsAtEnd(); ++bit){
    if (bit.Get()) nbPtsInBand++;
  }

  std::cout << "points in band: " << nbPtsInBand << std::endl;

  BandPointSampleVector& bandPointSampleVector = m_C[collision_time].bandPointSampleVector;
  bandPointSampleVector.resize(nbPtsInBand);

  unsigned int idx = 0;
  for( bit.GoToBegin(), lit.GoToBegin(); !bit.IsAtEnd(); ++bit, ++lit) {
    if (bit.Get()) {
      bandPointSampleVector[idx].index = bit.GetIndex();
      bandPointSampleVector[idx].label = lit.Get();

      MeasurementVectorType queryPoint;
      for(unsigned int ii=0; ii<InputDimension; ii++) queryPoint[ii] = bandPointSampleVector[idx].index[ii];
      for ( unsigned int i = 0; i < InputDimension; ++i ) origin[i] = queryPoint[i];
      distanceMetric->SetOrigin( origin );

      // querey for each kd-tree by label
      const float band_radius_thres = m_InterTrajectoryDistThres;
      bandPointSampleVector[idx].neighbor_labels.clear();
      bandPointSampleVector[idx].neighbor_distances.clear();
      bandPointSampleVector[idx].neighbor_weights.clear();

      for(unsigned int k=0; k<nbAffine; k++) {
        unsigned int label = k + 1;
        // K-Neighbor search
          typename TreeType::InstanceIdentifierVectorType neighbors;
          unsigned int numberOfNeighbors = 1;
          treeListForEachAffine[k]->Search( queryPoint, numberOfNeighbors, neighbors );

          // add the nearest point of label "=label" if the distance is below threshold
          // for ( unsigned int i = 0 ; i < neighbors.size() ; ++i ) { // only loop for once should be executed
          {
            unsigned int i = 0;
            double distance = distanceMetric->Evaluate( treeListForEachAffine[k]->GetMeasurementVector( neighbors[i] ) );
            // if (distance < band_radius_thres) {
            if (1) {
              bandPointSampleVector[idx].neighbor_labels.push_back(label);

              // point distance is negative when computing from points in the edge of the same label in the voronoi
              // point distance is positive when computing from points in the edge of the DIFFERENT label in the voronoi
              double polarized_distance = (label == bandPointSampleVector[idx].label) ? (-1*distance) : distance;
              bandPointSampleVector[idx].neighbor_distances.push_back(polarized_distance);

              float weight = compute_weight_from_distance_to_boundary(polarized_distance, m_InterTrajectoryDistThres);
              bandPointSampleVector[idx].neighbor_weights.push_back(weight);
              std::cout << "label=" << k+1 << " from " << queryPoint << " found " << treeListForEachAffine[k]->GetMeasurementVector( neighbors[i] ) << " d=" << polarized_distance << " w=" << weight << std::endl;
            }
          }
      }
      // normalize the weight
      unsigned int nbWeights = bandPointSampleVector[idx].neighbor_weights.size();
      if (nbWeights==0) std::cout << "Weired!! nbWeights = 0 for index =" << bandPointSampleVector[idx].index << " and label=" <<  bandPointSampleVector[idx].label  << std::endl;
      float wsum = 0.0;
      for(unsigned int i=0; i<nbWeights; i++) wsum += bandPointSampleVector[idx].neighbor_weights[i];
      std::cout << "Wsum = " << wsum << " for index =" << bandPointSampleVector[idx].index << " and label=" <<  bandPointSampleVector[idx].label  << std::endl;
      if (wsum <0.5) std::cout << "=======>Weired!! wsum <0.5 for index =" << bandPointSampleVector[idx].index << " and label=" <<  bandPointSampleVector[idx].label  << std::endl;
      for(unsigned int i=0; i<nbWeights; i++) bandPointSampleVector[idx].neighbor_weights[i] /= wsum;

      idx++;
    }
  }

}


template <class TScalar, unsigned int NDimensions>
void GSPolyAffineTransform<TScalar, NDimensions>::
GetTrajectoryVoronoiBand(const LabelImagePointerType &voronoi, unsigned int voronoi_band_radius, BinaryImagePointerType &band, BinaryImagePointerType &edge)
{
  // 1. get edge map of voronoi using sobel + threshold : dL
  // 2. dilate dL with radius = radius
  // 4. convert dL to list of each label
  // . for each point in dL, find the labels of all overlapping voronoi dilated band and store them using a list: for sparsity

  typedef SobelEdgeDetectionImageFilter<LabelImageType, FloatImageType>
  SobelEdgeDetectionImageFilterType;
  typename SobelEdgeDetectionImageFilterType::Pointer sobelFilter
  = SobelEdgeDetectionImageFilterType::New();
  sobelFilter->SetInput( voronoi );

  typedef itk::BinaryThresholdImageFilter<FloatImageType, BinaryImageType>
  BinaryThresholdImageFilterType;

  typename BinaryThresholdImageFilterType::Pointer thresholdFilter
  = BinaryThresholdImageFilterType::New();
  thresholdFilter->SetInput( sobelFilter->GetOutput() );
  thresholdFilter->SetLowerThreshold( 0.5 );
  thresholdFilter->SetInsideValue( 1 );
  thresholdFilter->SetOutsideValue( 0 );

  typedef itk::BinaryBallStructuringElement<typename BinaryImageType::PixelType, InputDimension> StructuringElementType;
  StructuringElementType structuringElement;
  std::cout << "voronoi_band_radius=" << voronoi_band_radius << std::endl;
  structuringElement.SetRadius( voronoi_band_radius );
  structuringElement.CreateStructuringElement();

  typedef itk::BinaryDilateImageFilter<BinaryImageType, BinaryImageType, StructuringElementType> BinaryDilateImageFilterType;

  typename BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  dilateFilter->SetInput( thresholdFilter->GetOutput() );
  dilateFilter->SetKernel( structuringElement );
  dilateFilter->SetDilateValue( 1 );
  dilateFilter->Update();


  band = dilateFilter->GetOutput();
  edge = thresholdFilter->GetOutput();

  PicslImageHelper::WriteImage<BinaryImageType>(band, "band.nii.gz");
  PicslImageHelper::WriteImage<BinaryImageType>(edge, "edge.nii.gz");

}

template <class TScalar, unsigned int NDimensions>
void GSPolyAffineTransform<TScalar, NDimensions>::ScatterLabeledPointsIntoTrajectoryImage(unsigned int collision_time, LabelImagePointerType &traj, SNIteratorType &it)
{
  unsigned int ind_t1 = m_CollisionTimeList[collision_time];
  unsigned int ind_t2 = m_CollisionTimeList[collision_time+1]-1;

  std::cout << m_CollisionTimeList[collision_time] << " --- " << m_CollisionTimeList[collision_time+1]-1 << " dkdjfkd " << m_ControlPointListTrajectory[1].size() <<  std::endl;

  unsigned int nbPts = m_ControlPointList.size();

  typename LabelImageType::RegionType region = traj->GetBufferedRegion();

  for(unsigned int i = 0; i < nbPts; i++) {
    unsigned int label = m_ControlPointLabelList[i];
    for(unsigned int t = ind_t1; t<= ind_t2; t++) {
      typename LabelImageType::IndexType index;
      typename LabelImageType::PointType pt;
      pt = m_ControlPointListTrajectory[i][t];
      m_LabelImage->TransformPhysicalPointToIndex(pt, index);
      it.SetLocation(index);

     typename SNIteratorType::Iterator ci;
      for (ci = it.Begin(); ! ci.IsAtEnd(); ++ci) {
        typename LabelImageType::IndexType innerIndex;
        innerIndex = index + ci.GetNeighborhoodOffset();
        if (region.IsInside(innerIndex))
          ci.Set(label);
      }
    }
  }
}


template <class TScalar, unsigned int NDimensions>
void GSPolyAffineTransform<TScalar, NDimensions>::PrecomputeTrajectoryAndCheckCollision2()
{
  // pre-compute all
  unsigned int nbAffine = m_LocalAffineList.size();
  unsigned int nb_T = m_NumberOfTimeSteps;
  unsigned int nbPts = m_ControlPointList.size();
  double inter_trajectory_dist2_thres = m_InterTrajectoryDistThres * m_InterTrajectoryDistThres;
  int cnt = 0;

  m_CollisionTimeList.clear();

  for(unsigned int j = 0; j < nbAffine; j++){
    m_LocalAffineList[j]->SetTimeStampMax(nb_T);
  }

  // initialize points trajectory [pts_idx][time] = single_point
  m_ControlPointListTrajectory.clear();
  m_ControlPointListTrajectory.resize(nbPts);
  for(unsigned int i=0; i < nbPts; i++) {
    m_ControlPointListTrajectory[i].resize(nb_T);
  }

  std::cout << "nbPts=" << nbPts << std::endl;
  std::cout << "nb_T=" << nb_T << std::endl;
  std::cout << "nbAffine=" << nbAffine << std::endl;

  std::vector<typename LocalAffineTransformType::AffineTransformPointer> cache_partial_transform;

  // non-collision time period is from m_CollisionTimeList[i] ~ m_CollisionTimeList[i+1]-1
  m_CollisionTimeList.push_back(0);
  // start track collision by time
  LabelImagePointerType collision_mask = LabelImageType::New();
  collision_mask->CopyInformation( m_LabelImage );
  collision_mask->SetRequestedRegion( m_LabelImage->GetRequestedRegion() );
  collision_mask->SetBufferedRegion( m_LabelImage->GetBufferedRegion() );
  collision_mask->Allocate();
  collision_mask->FillBuffer(0);
  for(unsigned int t = 0; t < nb_T; t++) {
    cache_partial_transform.clear();
    cache_partial_transform.resize(nbAffine);
    for(unsigned int j = 0; j < nbAffine; j++){
      cache_partial_transform[j] = m_LocalAffineList[j]->GetPartialTransform(0, t);
    }

    // extend trajectory to time=t for each point
    for(unsigned int i = 0; i < m_ControlPointList.size(); i++) {
      InputPointType pt = m_ControlPointList[i];
      LabelType j = m_ControlPointLabelList[i] - 1; // partial transform stored from 0 for mask label = 1
      InputPointType new_pt = cache_partial_transform[j]->TransformPoint(pt);
      m_ControlPointListTrajectory[i][t] = new_pt;
    }

    // computing pairwise distances seems very slow, try to use an image mask to hold all these points and fit a sphere with the radius

    unsigned int last_time = m_CollisionTimeList[m_CollisionTimeList.size()-1];
    bool b_collision = false;

    // Add New points in m_ControlPointListTrajectory[*][t] into the mask collision_mask
    for(unsigned int k=0; k<nbAffine && !b_collision; k++) {
      for(unsigned int i=0; i<nbPts && !b_collision; i++) {
        if (m_ControlPointLabelList[i] != k+1) continue; // label is k+1 not k
        cnt++;
        b_collision = !(AddSphereToLabelMask(collision_mask, m_ControlPointListTrajectory[i][t], m_InterTrajectoryDistThres, inter_trajectory_dist2_thres, k+1));
      }
    }

    if (b_collision && t>0) { // no need to add t=0 again, t=0 and collision simply means masks are too close from beginning!!
      collision_mask->FillBuffer(0);    // clear collsion_mask
      m_CollisionTimeList.push_back(t); // update collsion time
    }
  }

  m_CollisionTimeList.push_back(nb_T);

 // PicslImageHelper::WriteResultImage<LabelImageType>(collision_mask, "collsion_mask.nii.gz");
  std::cout << "add " << cnt << " spheres" << std::endl;
  std::cout << "collision time:" << std::endl;
  std::copy( m_CollisionTimeList.begin(), m_CollisionTimeList.end(), std::ostream_iterator<float>(std::cout, " "));
  std::cout <<  std::endl;
}


template <class TScalar, unsigned int NDimensions>
void GSPolyAffineTransform<TScalar, NDimensions>::FillBallShapedNeighborhoodIterator(const unsigned int radius, SNIteratorType &it){

  // Set active flags of iterator by looping an NeighborhoodIterator
  typedef BinaryBallStructuringElement<LabelType, InputDimension> StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();
  for(unsigned int i=0; i<structuringElement.Size(); i++) {
    if (structuringElement[i] > 0) {
      it.ActivateOffset(structuringElement.GetOffset(i));
    }
  }
  std::cout << "Now there are " << it.GetActiveIndexListSize() << " active indices." << std::endl;
  typename SNIteratorType::IndexListType indexList = it.GetActiveIndexList();
  typename SNIteratorType::IndexListType::const_iterator listIterator = indexList.begin();
  while (listIterator != indexList.end())
   {
   std::cout << *listIterator << " ";
   ++listIterator;
   }
  std::cout << std::endl;

  std::cout << "ball neighborhood:" << structuringElement << std::endl;


  // return it;
}


template <class TScalar, unsigned int NDimensions>
bool GSPolyAffineTransform<TScalar, NDimensions>::AddSphereToLabelMask(LabelImagePointerType &collision_image, const InputPointType &pt, const unsigned int radius, const double inter_trajectory_dist2_thres, const unsigned int labelK){

  typename LabelImageType::RegionType bbox;
  typename LabelImageType::IndexType center, corner;
  typename itk::ContinuousIndex<double, InputDimension> ccenter;
  typename LabelImageType::OffsetType offset;
  typename LabelImageType::SizeType size;

  offset.Fill(radius+1);
  collision_image->TransformPhysicalPointToIndex(pt, center);
  collision_image->TransformPhysicalPointToContinuousIndex(pt, ccenter);
  corner = center - offset;
  size.Fill(2*radius + 1);
  bbox.SetIndex(corner);
  bbox.SetSize(size);

  bool bcrop = bbox.Crop(collision_image->GetBufferedRegion());

//  std::cout << " bbox=" << bbox << " bcrop=" << bcrop << std::endl;

  if (bcrop) {
    itk::ImageRegionIteratorWithIndex<LabelImageType> iter(collision_image, bbox);
    for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
      typename LabelImageType::IndexType position;
      position = iter.GetIndex();
      double distanceSquared = 0.0;
      for ( unsigned int i = 0; i < InputDimension; i++ )
        distanceSquared += (position[i] -ccenter[i]) * (position[i] -ccenter[i]);


      //    std::cout << "distanceSquared="<<distanceSquared << " position="<<position << " ccenter=" << ccenter << " radius=" << radius
      //        << " inter_trajectory_dist2_thres=" << inter_trajectory_dist2_thres
      //        << "exist label=" << iter.Get() << " k=" << labelK <<  std::endl;

      if ( distanceSquared <= inter_trajectory_dist2_thres ) {
        unsigned int existing_label = iter.Get();
        iter.Set(labelK);
        if ( existing_label > 0 && existing_label != labelK) // already existed another label value
          {
          return false;
          }
      }
      //    std::cout << "new label=" << iter.Get() << " k=" << labelK <<  std::endl;
    }
  }
  return true;
}



template <class TScalar, unsigned int NDimensions>
void GSPolyAffineTransform<TScalar, NDimensions>::PrecomputeTrajectoryAndCheckCollision()
{
  // pre-compute all
  unsigned int nbAffine = m_LocalAffineList.size();
  unsigned int nb_T = 50;
  unsigned int nbPts = m_ControlPointList.size();

  m_CollisionTimeList.clear();

  for(unsigned int j = 0; j < nbAffine; j++){
    m_LocalAffineList[j]->SetTimeStampMax(nb_T);
  }

  // initialize points trajectory [pts_idx][time] = single_point
  m_ControlPointListTrajectory.clear();
  m_ControlPointListTrajectory.resize(nbPts);
  for(unsigned int i=0; i < nbPts; i++) {
    m_ControlPointListTrajectory[i].resize(nb_T);
  }

  std::cout << "nbPts=" << nbPts << std::endl;
  std::cout << "nb_T=" << nb_T << std::endl;
  std::cout << "nbAffine=" << nbAffine << std::endl;

  std::vector<typename LocalAffineTransformType::AffineTransformPointer> cache_partial_transform;

  // non-collision time period is from m_CollisionTimeList[i] ~ m_CollisionTimeList[i+1]-1
  m_CollisionTimeList.push_back(0);
  // start track collision by time

  for(unsigned int t = 0; t < nb_T; t++) {
    cache_partial_transform.clear();
    cache_partial_transform.resize(nbAffine);
    for(unsigned int j = 0; j < nbAffine; j++){
      cache_partial_transform[j] = m_LocalAffineList[j]->GetPartialTransform(0, t);
    }

    // extend trajectory to time=t for each point
    for(unsigned int i = 0; i < m_ControlPointList.size(); i++) {
      InputPointType pt = m_ControlPointList[i];
      LabelType j = m_ControlPointLabelList[i] - 1; // partial transform stored from 0 for mask label = 1
      InputPointType new_pt = cache_partial_transform[j]->TransformPoint(pt);
      m_ControlPointListTrajectory[i][t] = new_pt;
    }

    // check from the last collision time point
    unsigned int last_time = m_CollisionTimeList[m_CollisionTimeList.size()-1];

    double d2;
    double inter_trajectory_dist2_thres = m_InterTrajectoryDistThres * m_InterTrajectoryDistThres;
    bool b_collision = false;

    if (t == last_time) {
      // start computing the first point in collision at last_time vs last time;
      for(unsigned int i=0; i<nbPts && !b_collision; i++) {
        for(unsigned int j=0; j<nbPts && !b_collision; j++) {
          if (m_ControlPointLabelList[i] == m_ControlPointLabelList[j])
            continue;
          d2 = m_ControlPointListTrajectory[i][t].SquaredEuclideanDistanceTo(m_ControlPointListTrajectory[j][t]);
          if (d2 < inter_trajectory_dist2_thres) {
            b_collision = true;
          }
        }
      }
    }

    if (t > last_time) {
      for(unsigned int i=0; i<nbPts && !b_collision; i++) {
        for(unsigned int j=0; j<nbPts && !b_collision; j++) {
          if (m_ControlPointLabelList[i] == m_ControlPointLabelList[j])
            continue;

          unsigned int ti, tj;
          // compute distance of trajectory incrementally by time
          // t vs t, t vs [last_time, t-1], [last_time, t-1] vs t,

          tj = t;
          for(ti=last_time; ti<t && !b_collision; ti++) {
            d2 = m_ControlPointListTrajectory[i][ti].SquaredEuclideanDistanceTo(m_ControlPointListTrajectory[j][tj]);
            if (d2 < inter_trajectory_dist2_thres) {
              b_collision = true;
            }
          }

          ti = t;
          for(tj=last_time; tj<t && !b_collision; tj++) {
            d2 = m_ControlPointListTrajectory[i][ti].SquaredEuclideanDistanceTo(m_ControlPointListTrajectory[j][tj]);
            if (d2 < inter_trajectory_dist2_thres) {
              b_collision = true;
            }
          }

          ti = t; tj = t;
          d2 = m_ControlPointListTrajectory[i][ti].SquaredEuclideanDistanceTo(m_ControlPointListTrajectory[j][tj]);
          if (d2 < inter_trajectory_dist2_thres) {
            b_collision = true;
          }

        }
      }
    }

    if (b_collision) {
      m_CollisionTimeList.push_back(t);
    }
  }

  m_CollisionTimeList.push_back(nb_T+1);

  std::cout << "collision time:" << std::endl;
  std::copy( m_CollisionTimeList.begin(), m_CollisionTimeList.end(), std::ostream_iterator<float>(std::cout, " "));
}


template <class TScalar, unsigned int NDimensions>
void GSPolyAffineTransform<TScalar, NDimensions>::ConvertAffineTransformsToLogarithm()
{
	m_LocalAffineList.clear();
	unsigned int nbAffine = this->GetNumberOfTransforms();
	m_LocalAffineList.resize(nbAffine);

	for(unsigned i=0; i<nbAffine; i++) {
		typename LocalAffineTransformType::Pointer laff = LocalAffineTransformType::New();
		typename MatrixOffsetTransformBaseType::Pointer aff = dynamic_cast<MatrixOffsetTransformBaseType *> (this->GetNthTransform(i).GetPointer());
		
		//laff should be centerd at zero
		laff->SetMatrix(aff->GetMatrix());
		laff->SetOffset(aff->GetOffset());

		// get velocity matrix
		laff->GetVelocityMatrix();

		std::cout << "[" << i << "]" << aff << std::endl;
		std::cout << "==>[" << i << "]" << laff << std::endl;

		m_LocalAffineList[i] = laff;

	}



}



template <class TScalar, unsigned int NDimensions>
void GSPolyAffineTransform<TScalar, NDimensions>::  SampleFromLabelMask(const LabelImagePointerType &labelImage, const SampleSpacingType &sample_spacing, PointVectorType &controlPointList, PointLabelVectorType &controlPointLabelList)
{
	typedef typename LabelImagePointerType::ObjectType LabelImageType;

	typename LabelImageType::RegionType region = labelImage->GetLargestPossibleRegion();
	typename LabelImageType::SizeType imgsz = region.GetSize();
	
	// loop over x/y/z, only support 2D/3D
	unsigned int dimz = (InputDimension > 2) ? imgsz[2] : 1;
	unsigned int dimy = imgsz[1];
	unsigned int dimx = imgsz[0];

	std::vector<unsigned int> spacing;
	spacing.resize(3);
	for(unsigned int i=0; i<3; i++) 
		spacing[i] = 1;
	for(unsigned int i=0; i<InputDimension; i++)
		spacing[i] = sample_spacing[i];

	unsigned int nb_max_points = ceil(dimz * dimy * dimx / (sample_spacing[0] * sample_spacing[1] * sample_spacing[2]));

	controlPointList.reserve(nb_max_points);
	controlPointLabelList.reserve(nb_max_points);

	typename LabelImageType::IndexType index;
	for(unsigned int z = 0; z < dimz; z+= spacing[2]) {
		if (InputDimension > 2)	index[2] = z;
		for(unsigned int y = 0; y < dimy; y+= spacing[1]) {
			index[1] = y;
			for(unsigned int x = 0; x < dimx; x+= spacing[0]) {
				index[0] = x;
				LabelType label = labelImage->GetPixel(index);

				std::cout << "[x,y,z]" << x << "," << y << ',' << z << "=" << label << " dimz=" << dimz << std::endl;

				if ( label > 0) {
					InputPointType point;
					labelImage->TransformIndexToPhysicalPoint(index, point);
					controlPointList.push_back(point);
					controlPointLabelList.push_back(label);
				}
			}
		}
	}
}



}


#endif
