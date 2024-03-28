set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_SOURCE_DIR}/core/mpdBase
  ${CMAKE_SOURCE_DIR}/core/mpdPid
  ${CMAKE_SOURCE_DIR}/detectors/zdc
  ${CMAKE_SOURCE_DIR}/physics
  ${CMAKE_SOURCE_DIR}/physics/evPlane
  ${CMAKE_SOURCE_DIR}/physics/evFlowEP
  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/kalman
  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/lheTrack
  ${CMAKE_SOURCE_DIR}/simulation/mcStack
)

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})
include_Directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
  ${ROOT_LIBRARY_DIR}
  ${FAIRROOT_LIBRARY_DIR}
)

link_directories(${LINK_DIRECTORIES})

set(SRCS
  MpdFlowEventPlane.cxx
  MpdFlowEventPlaneParams.cxx
)

Set(HEADERS)
Set(LINKDEF MpdFlowEventPlaneLinkDef.h)
Set(LIBRARY_NAME MpdFlowEventPlane)
Set(DEPENDENCIES)

GENERATE_LIBRARY()
