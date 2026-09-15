set ( qcbiasfilter_files
  QcBiasFilterBase.h
  QcBiasFilterBase.cc
  CloudDiffFilter.h
  CloudDiffFilter.cc
)

PREPEND( _p_qcbiasfilter_files     "qcbiasfilters"            ${qcbiasfilter_files} )

set ( qcbiasfilter_src_files
      ${_p_qcbiasfilter_files}
      PARENT_SCOPE
)
