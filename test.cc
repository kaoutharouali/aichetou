/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */
/*
 *   Copyright (c) 2020 Centre Tecnologic de Telecomunicacions de Catalunya (CTTC)
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2 as
 *   published by the Free Software Foundation;
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/**
 * \ingroup examples
 * \file cttc-nr-v2x-demo-simple.cc
 * \brief A cozy, simple, NR V2X demo (in a tutorial style)
 *
 * This example describes how to setup an NR sidelink out-of-coverage simulation
 * using the 3GPP channel model from TR 37.885. This example simulates a
 * simple topology consist of 2 UEs, where UE-1 transmits, and UE-2 receives.
 *
 * Have a look at the possible parameters to know what you can configure
 * through the command line.
 *
 * With the default configuration, the example will create a flow that will
 * go through a subband or bandwidth part. For that,
 * specifically, one band with a single CC, and one bandwidth part is used.
 *
 * The example will print on-screen the average Packet Inter-Reception (PIR)
 * type 2 computed as defined in 37.885. Moreover, it saves MAC and PHY layer
 * traces in a sqlite3 database using ns-3 stats module. Moreover, since there
 * is only one transmitter in the scenario, sensing is by default not enabled.
 *
 * \code{.unparsed}
$ ./waf --run "cttc-nr-v2x-demo-simple --help"
    \endcode
 *
 */

/*
 * Include part. Often, you will have to include the headers for an entire module;
 * do that by including the name of the module you need with the suffix "-module.h".
 */
#include "ns3/core-module.h"
#include "ns3/config-store.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/applications-module.h"
#include "ns3/mobility-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/nr-module.h"
#include "ns3/lte-module.h"
#include "ns3/stats-module.h"
#include "ns3/config-store-module.h"
#include "ns3/log.h"
#include "ue-mac-pscch-tx-output-stats.h"
#include "ue-mac-pssch-tx-output-stats.h"
#include "ue-phy-pscch-rx-output-stats.h"
#include "ue-phy-pssch-rx-output-stats.h"
#include "ue-to-ue-pkt-txrx-output-stats.h"
#include "ns3/olsr-module.h"
#include "ns3/aodv-module.h"
#include "ns3/antenna-module.h"
#include "ns3/olsr-helper.h"
#include "ns3/olsr-routing-protocol.h"
#include "ns3/ipv4-static-routing-helper.h"
#include "ns3/ipv4-list-routing-helper.h"
#include "v2x-kpi.h"
#include "ue-rlc-rx-output-stats.h"
#include "ns3/flow-monitor-module.h"
#include <iomanip>


/*
 * Use, always, the namespace ns3. All the NR classes are inside such namespace.
 */
using namespace ns3;

/*
 * With this line, we will be able to see the logs of the file by enabling the
 * component "CttcNrV2xDemoSimple", in this way:
 *
 * $ export NS_LOG="CttcNrV2xDemoSimple=level_info|prefix_func|prefix_time"
 */
NS_LOG_COMPONENT_DEFINE ("test");

/*
 * Global variables to count TX/RX packets and bytes.
 */

 void NotifySlPscchScheduling (UeMacPscchTxOutputStats *pscchStats, const SlPscchUeMacStatParameters pscchStatsParams)
 {
   pscchStats->Save (pscchStatsParams);
 }

 void NotifySlPsschScheduling (UeMacPsschTxOutputStats *psschStats, const SlPsschUeMacStatParameters psschStatsParams)
 {
   psschStats->Save (psschStatsParams);
 }
 void NotifySlPscchRx (UePhyPscchRxOutputStats *pscchStats, const SlRxCtrlPacketTraceParams pscchStatsParams)
 {
   pscchStats->Save (pscchStatsParams);
 }


 void NotifySlPsschRx (UePhyPsschRxOutputStats *psschStats, const SlRxDataPacketTraceParams psschStatsParams)
 {
   psschStats->Save (psschStatsParams);
 }

 void NotifySlRlcPduRx (UeRlcRxOutputStats *stats, uint64_t imsi, uint16_t rnti, uint16_t txRnti, uint8_t lcid, uint32_t rxPduSize, double delay)
 {
   stats->Save (imsi, rnti, txRnti, lcid, rxPduSize, delay);
 }

 void
 UePacketTraceDb (UeToUePktTxRxOutputStats *stats, Ptr<Node> node, const Address &localAddrs,
                  std::string txRx, Ptr<const Packet> p, const Address &srcAddrs,
                  const Address &dstAddrs, const SeqTsSizeHeader &seqTsSizeHeader)
 {
   uint32_t nodeId = node->GetId ();
   uint64_t imsi = node->GetDevice (0)->GetObject<NrUeNetDevice> ()->GetImsi ();
   uint32_t seq = seqTsSizeHeader.GetSeq ();
   uint32_t pktSize = p->GetSize () + seqTsSizeHeader.GetSerializedSize ();

   stats->Save (txRx, localAddrs, nodeId, imsi, pktSize, srcAddrs, dstAddrs, seq);
 }


 /**
  * \brief Install highway mobility
  * \param totalLanes Total number of lanes
  * \param numVehiclesPerLane number of vehicles per lane (only odd numbers)
  * \param interVehicleDist The distance between the vehicles in a same lane
  * \param interLaneDist The distance between the lanes, i.e., the distance
  *        between the vehicles of two different lanes
  * \param speed The speed of the vehicles
  * \return A node container containing the vehicle UEs with their mobility model
  *         installed
  */
/* NodeContainer
 //InstallHighwayMobility (uint16_t totalLanes, uint16_t numVehiclesPerLane, uint16_t interVehicleDist, uint16_t interLaneDist, double speed)
 //{
  // NS_ABORT_MSG_IF ((numVehiclesPerLane * totalLanes) % totalLanes != 0, "All the lanes must have equal number of UEs");

   NodeContainer ueNodes;

   ueNodes.Create (numVehiclesPerLane * totalLanes);

   std::cout << "Total UEs = " << ueNodes.GetN () << std::endl;

   Ptr<GridPositionAllocator>  gridPositionAllocator;
   gridPositionAllocator = CreateObject<GridPositionAllocator> ();
   gridPositionAllocator->SetAttribute ("GridWidth", UintegerValue (numVehiclesPerLane));
   gridPositionAllocator->SetAttribute ("MinX", DoubleValue (0.0));
   gridPositionAllocator->SetAttribute ("MinY", DoubleValue (0.0));
   gridPositionAllocator->SetAttribute ("Z", DoubleValue (1.6));
   gridPositionAllocator->SetAttribute ("DeltaX", DoubleValue (interVehicleDist));
   gridPositionAllocator->SetAttribute ("DeltaY", DoubleValue (interLaneDist));
   gridPositionAllocator->SetAttribute ("LayoutType", EnumValue (GridPositionAllocator::ROW_FIRST));

   MobilityHelper mobility;
   mobility.SetMobilityModel ("ns3::ConstantVelocityMobilityModel");
   mobility.SetPositionAllocator (gridPositionAllocator);

   mobility.Install (ueNodes);

   for (int i = 0; i < totalLanes * numVehiclesPerLane; i++)
     {
       ueNodes.Get (i)->GetObject<ConstantVelocityMobilityModel> ()->SetVelocity (Vector (speed, 0.0, 0.0));
     }

   return ueNodes;
// }
*/
 /**
  * \brief Write the gnuplot script to plot the initial positions of the vehicle UEs
  * \param posFilename The name of the file, which this script needs to read to plot positions
  */
 /*void
 WriteInitPosGnuScript (std::string posFilename)
 {
   std::ofstream outFile;
   std::string filename = "gnu-script-" + posFilename;
   outFile.open (filename.c_str (), std::ios_base::out | std::ios_base::trunc);
   if (!outFile.is_open ())
     {
       NS_LOG_ERROR ("Can't open file " << filename);
       return;
     }

   std::string pngFileName;
   pngFileName = posFilename.substr (0, posFilename.rfind ("."));
   outFile << "set terminal png" << std::endl;
   outFile << "set output \"" << pngFileName << ".png\"" << std::endl;
   outFile << "set style line 1 lc rgb 'black' ps 2 pt 7" << std::endl;
   outFile << "unset key" << std::endl;
   outFile << "set grid" << std::endl;
   outFile << "plot \"" << posFilename << "\" using 3:4 with points ls 1";
   outFile.close ();
 }
*/
 /**
  * \brief Print vehicle UEs initial position to a file
  * \param filename Name of the file in which to write the positions
  */
/* void
 PrintUeInitPosToFile (std::string filename)
 {
   std::ofstream outFile;
   outFile.open (filename.c_str (), std::ios_base::out | std::ios_base::trunc);
   if (!outFile.is_open ())
     {
       NS_LOG_ERROR ("Can't open file " << filename);
       return;
     }
   for (NodeList::Iterator it = NodeList::Begin (); it != NodeList::End (); ++it)
     {
       Ptr<Node> node = *it;
       int nDevs = node->GetNDevices ();
       for (int j = 0; j < nDevs; j++)
         {
           Ptr<NrUeNetDevice> uedev = node->GetDevice (j)->GetObject <NrUeNetDevice> ();
           if (uedev)
             {
               Vector pos = node->GetObject<MobilityModel> ()->GetPosition ();
               outFile << node->GetId () << " " << uedev->GetImsi () << " " << pos.x << " " << pos.y << std::endl;
             }
         }
     }

   WriteInitPosGnuScript (filename);
 }
*/
 /**
  * \brief Record mobility of the vehicle UEs every second
  * \param FirstWrite If this flag is true, write from scratch, otherwise, append to the file
  * \param fileName Name of the file in which to write the positions
  */
/* void RecordMobility (bool FirstWrite, std::string fileName)
 {
   std::ofstream outFile;
   if (FirstWrite == true)
     {
       outFile.open (fileName.c_str (), std::ios_base::out);
       FirstWrite = false;
     }
   else
     {
       outFile.open (fileName.c_str (), std::ios_base::app);
       outFile << std::endl;
       outFile << std::endl;
     }

   for (NodeList::Iterator it = NodeList::Begin (); it != NodeList::End (); ++it)
     {
       Ptr<Node> node = *it;
       int nDevs = node->GetNDevices ();
       for (int j = 0; j < nDevs; j++)
         {
           Ptr<NrUeNetDevice> uedev = node->GetDevice (j)->GetObject <NrUeNetDevice> ();
           if (uedev)
             {
               Vector pos = node->GetObject<MobilityModel> ()->GetPosition ();
               outFile << Simulator::Now ().GetSeconds () << " " << node->GetId () << " " << uedev->GetImsi () << " " << pos.x << " " << pos.y << std::endl;
             }
         }
     }

   Simulator::Schedule (Seconds (1), &RecordMobility, FirstWrite, fileName);
 }
*/
 /**
  * \brief Write a gnuplot script to generate gif of the vehicle UEs mobility
  * \param MobilityFileName The name of the file, which this script needs to read to plot positions
  * \param simTime The simulation time
  * \param speed The speed of the vehicles
  * \param firstUeNode The node pointer to the first UE in the simulation
  * \param lastUeNode The node pointer to the last UE in the simulation
  */
/* void WriteGifGnuScript (std::string MobilityFileName, Time simTime, double speed, Ptr<Node> firstUeNode, Ptr<Node> lastUeNode)
 {
   std::ofstream outFile;
   std::string fileName = "gif-script-" + MobilityFileName;
   outFile.open (fileName.c_str (), std::ios_base::out | std::ios_base::trunc);
   if (!outFile.is_open ())
     {
       NS_LOG_ERROR ("Can't open file " << fileName);
       return;
     }
   outFile << "set term gif animate delay 100" << std::endl;
   std::string gifFileName;
   gifFileName = MobilityFileName.substr (0, MobilityFileName.rfind ("."));
   outFile << "set output \"" << gifFileName << ".gif" << "\"" << std::endl;
   outFile << "unset key" << std::endl;
   outFile << "set grid" << std::endl;

   Vector firstNodePos = firstUeNode->GetObject<MobilityModel> ()->GetPosition ();
   Vector LastNodePos = lastUeNode->GetObject<MobilityModel> ()->GetPosition ();
   double xRangeLower = firstNodePos.x - 10.0;
   double xRangeUpper = simTime.GetSeconds () * speed + LastNodePos.x;
   double yRangeLower = firstNodePos.y - 10.0;
   double yRangeUpper = LastNodePos.y + 10.0;
   outFile << "set xrange [" << xRangeLower << ":" << xRangeUpper << "]" << std::endl;
   outFile << "set yrange [" << yRangeLower << ":" << yRangeUpper << "]" << std::endl;
   outFile << "do for [i=0:" << simTime.GetSeconds () - 1 << "] {plot \"" << MobilityFileName << "\" index i using 4:5}" << std::endl;
 }
*/
 /**
  * \brief Get sidelink bitmap from string
  * \param slBitMapString The sidelink bitmap string
  * \param slBitMapVector The vector passed to store the converted sidelink bitmap
  */
 /*void
 GetSlBitmapFromString (std::string slBitMapString, std::vector <std::bitset<1> > &slBitMapVector)
 {
   static std::unordered_map<std::string, uint8_t> lookupTable =
   {
     { "0", 0 },
     { "1", 1 },
   };

   std::stringstream ss (slBitMapString);
   std::string token;
   std::vector<std::string> extracted;

   while (std::getline (ss, token, '|'))
     {
       extracted.push_back (token);
     }

   for (const auto & v : extracted)
     {
       if (lookupTable.find (v) == lookupTable.end ())
         {
           NS_FATAL_ERROR ("Bit type " << v << " not valid. Valid values are: 0 and 1");
         }
       slBitMapVector.push_back (lookupTable[v] & 0x01);
     }
 }*/
/**
 * \brief Method to listen the packet sink application trace Rx.
 * \param packet The packet
 * \param from The address of the transmitting node
 */
//void ReceivePacket (Ptr<const Packet> packet, const Address & from)
//{
 //NS_UNUSED (from);
//rxByteCounter += packet->GetSize();
//  rxPktCounter++;
//}

/**
 * \brief Method to listen the transmitting application trace Tx.
 * \param packet The packet
 */

// change to code socket add

uint32_t txPktCounter = 0;
uint32_t txByteCounter = 0;

void TransmitPacket (Ptr<const Packet> packet)
{
  txByteCounter += packet->GetSize();
  txPktCounter++;
}

class RoutingStats
{
public:
  /**
   * \brief Constructor
   */
  RoutingStats ();

  /**
   * \brief Returns the number of bytes received
   * \return the number of bytes received
   */
  uint32_t GetRxBytes ();

  /**
   * \brief Returns the cumulative number of bytes received
   * \return the cumulative number of bytes received
   */
  uint32_t GetCumulativeRxBytes ();

  /**
   * \brief Returns the count of packets received
   * \return the count of packets received
   */
  uint32_t GetRxPkts ();

  /**
   * \brief Returns the cumulative count of packets received
   * \return the cumulative count of packets received
   */
  uint32_t GetCumulativeRxPkts ();

  /**
   * \brief Increments the number of (application-data)
   * bytes received, not including MAC/PHY overhead
   * \param rxBytes the number of bytes received
   */
  void IncRxBytes (uint32_t rxBytes);

  /**
   * \brief Increments the count of packets received
   */
  void IncRxPkts ();

  /**
   * \brief Sets the number of bytes received.
   * \param rxBytes the number of bytes received
   */
  void SetRxBytes (uint32_t rxBytes);

  /**
   * \brief Sets the number of packets received
   * \param rxPkts the number of packets received
   */
  void SetRxPkts (uint32_t rxPkts);

  /**
   * \brief Returns the number of bytes transmitted
   * \return the number of bytes transmitted
   */
  uint32_t GetTxBytes ();

  /**
   * \brief Returns the cumulative number of bytes transmitted
   * \return the cumulative number of bytes transmitted
   */
  uint32_t GetCumulativeTxBytes ();

  /**
   * \brief Returns the number of packets transmitted
   * \return the number of packets transmitted
   */
  uint32_t GetTxPkts ();

  /**
   * \brief Returns the cumulative number of packets transmitted
   * \return the cumulative number of packets transmitted
   */
  uint32_t GetCumulativeTxPkts ();

  /**
   * \brief Increment the number of bytes transmitted
   * \param txBytes the number of additional bytes transmitted
   */
  void IncTxBytes (uint32_t txBytes);

  /**
   * \brief Increment the count of packets transmitted
   */
  void IncTxPkts ();

  /**
   * \brief Sets the number of bytes transmitted
   * \param txBytes the number of bytes transmitted
   */
  void SetTxBytes (uint32_t txBytes);

  /**
   * \brief Sets the number of packets transmitted
   * \param txPkts the number of packets transmitted
   */
  void SetTxPkts (uint32_t txPkts);

private:
  uint32_t m_RxBytes; ///< reeive bytes
  uint32_t m_cumulativeRxBytes; ///< cumulative receive bytes
  uint32_t m_RxPkts; ///< receive packets
  uint32_t m_cumulativeRxPkts; ///< cumulative receive packets
  uint32_t m_TxBytes; ///< transmit bytes
  uint32_t m_cumulativeTxBytes; ///< cumulative transmit bytes
  uint32_t m_TxPkts; ///< transmit packets
  uint32_t m_cumulativeTxPkts; ///< cumulative transmit packets
};

RoutingStats::RoutingStats ()
  : m_RxBytes (0),
    m_cumulativeRxBytes (0),
    m_RxPkts (0),
    m_cumulativeRxPkts (0),
    m_TxBytes (0),
    m_cumulativeTxBytes (0),
    m_TxPkts (0),
    m_cumulativeTxPkts (0)
{
}

uint32_t
RoutingStats::GetRxBytes ()
{
  return m_RxBytes;
}

uint32_t
RoutingStats::GetCumulativeRxBytes ()
{
  return m_cumulativeRxBytes;
}

uint32_t
RoutingStats::GetRxPkts ()
{
  return m_RxPkts;
}

uint32_t
RoutingStats::GetCumulativeRxPkts ()
{
  return m_cumulativeRxPkts;
}

void
RoutingStats::IncRxBytes (uint32_t rxBytes)
{
  m_RxBytes += rxBytes;
  m_cumulativeRxBytes += rxBytes;
}

void
RoutingStats::IncRxPkts ()
{
  m_RxPkts++;
  m_cumulativeRxPkts++;
}

void
RoutingStats::SetRxBytes (uint32_t rxBytes)
{
  m_RxBytes = rxBytes;
}

void
RoutingStats::SetRxPkts (uint32_t rxPkts)
{
  m_RxPkts = rxPkts;
}

uint32_t
RoutingStats::GetTxBytes ()
{
  return m_TxBytes;
}

uint32_t
RoutingStats::GetCumulativeTxBytes ()
{
  return m_cumulativeTxBytes;
}

uint32_t
RoutingStats::GetTxPkts ()
{
  return m_TxPkts;
}

uint32_t
RoutingStats::GetCumulativeTxPkts ()
{
  return m_cumulativeTxPkts;
}

void
RoutingStats::IncTxBytes (uint32_t txBytes)
{
  m_TxBytes += txBytes;
  m_cumulativeTxBytes += txBytes;
}

void
RoutingStats::IncTxPkts ()
{
  m_TxPkts++;
  m_cumulativeTxPkts++;
}

void
RoutingStats::SetTxBytes (uint32_t txBytes)
{
  m_TxBytes = txBytes;
}

void
RoutingStats::SetTxPkts (uint32_t txPkts)
{
  m_TxPkts = txPkts;
}



class RoutingExperiment
{
public:
  RoutingExperiment ();
  void Run (int nSinks, double txPower, std::string CSVfileName);
  //static void SetMACParam (ns3::NetDeviceContainer & devices,
  //                                 int slotDistance);
  std::string CommandSetup (int argc, char **argv);

private:
  Ptr<Socket> SetupPacketReceive (Ipv4Address addr, Ptr<Node> node);
  void ReceivePacket (Ptr<Socket> socket);
  void CheckThroughput ();
  void ComputePir (Ptr<const Packet> packet, const Address &from);
  void OnOffTrace (std::string context, Ptr<const Packet> packet);
  //void TransmitPacket (Ptr<Packet> packet);

  RoutingStats & GetRoutingStats ();

  uint32_t port;
  uint32_t bytesTotal;
  uint32_t packetsReceived;
  uint32_t PktsSent;



  //uint32_t txPktCounter; //!< Global variable to count RX bytes
  //uint32_t txByteCounter; //!< Global variable to count TX bytes

  uint64_t pirCounter = 0; //!< counter to count how many time we computed the PIR. It is used to compute average PIR
  Time lastPktRxTime; //!< Global variable to store the RX time of a packet
  Time pir;

  std::string m_CSVfileName;
  int m_nSinks;
  std::string m_protocolName;
  double m_txp;
  //double pdr;
  bool m_traceMobility;
  uint32_t m_protocol;
  uint16_t m_numVehiclesPerLane;
  bool logging = false;

  RoutingStats routingStats; ///< routing statistics

   //bool enableOneTxPerLane = true;

   // uint16_t numVehiclesPerLane = 5;
  //uint16_t numLanes = 3;
   //double totalTime = 100;
   bool printRoutes;
  // Traffic parameters (that we will use inside this script:)
  bool useIPv6 = false; // default IPV4
  uint32_t udpPacketSizeBe = 300;
  double dataRateBe = 24; //16 kilobits per second

  //double speed = 38.88889; //meter per second, default 140 km/h

  // Simulation parameters.
  Time simTime = Seconds (10);
  //Sidelink bearers activation time
  Time slBearersActivationTime = Seconds (2.0);

  //flags to generate gnuplot plotting scripts
  bool generateInitialPosGnuScript = false;
  bool generateGifGnuScript = false;

  // NR parameters. We will take the input from the command line, and then we
  // will pass them inside the NR module.
  uint16_t numerologyBwpSl = 0;
  double centralFrequencyBandSl = 5.89e9; // band n47  TDD //Here band is analogous to channel
  uint16_t bandwidthBandSl = 100; //Multiple of 100 KHz; 400 = 40 MHz
  double txPower = 23; //dBm
  //double pdr = 0;
  double pDr = 0;
  // Where we will store the output files.
  std::string simTag = "default";
  std::string outputDir = "./";

};

RoutingExperiment::RoutingExperiment ()
  : port (9),
    bytesTotal (0),
     PktsSent (0),
    //pdr(0),
    //txPktCounter (0), //!< Global variable to count RX bytes
    //txByteCounter (0), //!< Global variable to count TX bytes
    packetsReceived (0),
    m_CSVfileName ("test-routing.output.csv"),
    m_traceMobility (false),
    //m_numVehiclesPerLane (nvLane),
    m_protocol (1) // olsr
{
}

/**
 * Print a received routing packet on a string
 * \param socket Rx socket
 * \param packet Rx packet
 * \param srcAddress source address
 * \return the built string
 */
static inline std::string
PrintReceivedPacket (Ptr<Socket> socket, Ptr<Packet> packet, Address senderAddress)
{
  std::ostringstream oss;

  oss << Simulator::Now ().GetSeconds () << " " << socket->GetNode ()->GetId ();

  if (InetSocketAddress::IsMatchingType (senderAddress))
    {
      InetSocketAddress addr = InetSocketAddress::ConvertFrom (senderAddress);
      oss << " received one packet from " << addr.GetIpv4 ();
    }
  else
    {
      oss << " received one packet!";
    }
  return oss.str ();
}
/*
void SavePositionPerIP (V2xKpi *v2xKpi)
{
  for (NodeList::Iterator it = NodeList::Begin (); it != NodeList::End (); ++it)
      {
        Ptr<Node> node = *it;
        int nDevs = node->GetNDevices ();
        for (int j = 0; j < nDevs; j++)
          {
            Ptr<NrUeNetDevice> uedev = node->GetDevice (j)->GetObject <NrUeNetDevice> ();
            if (uedev)
              {
                Ptr<Ipv4L3Protocol> ipv4Protocol =  node->GetObject<Ipv4L3Protocol> ();
                Ipv4InterfaceAddress addresses = ipv4Protocol->GetAddress (1,0);
                std::ostringstream ueIpv4Addr;
                ueIpv4Addr.str ("");
                ueIpv4Addr << addresses.GetLocal ();
                Vector pos = node->GetObject<MobilityModel> ()->GetPosition ();
                v2xKpi->FillPosPerIpMap (ueIpv4Addr.str (), pos);
              }
          }
      }
}
*/
void
RoutingExperiment::ReceivePacket (Ptr<Socket> socket)
{
  Ptr<Packet> packet;
  //txByteCounter += packet->GetSize();
  Address senderAddress;
  while ((packet = socket->RecvFrom (senderAddress)))
    {
      //bytesTotal += packet->GetSize ();
      uint32_t RxRoutingBytes = packet->GetSize ();
      GetRoutingStats ().IncRxBytes (RxRoutingBytes);
      GetRoutingStats ().IncRxPkts ();

    //  packetsReceived++;
      NS_LOG_UNCOND (PrintReceivedPacket (socket, packet, senderAddress));
    }
}

void
RoutingExperiment::OnOffTrace (std::string context, Ptr<const Packet> packet)
{
  uint32_t pktBytes = packet->GetSize ();
  routingStats.IncTxBytes (pktBytes);
  GetRoutingStats ().IncTxPkts ();
}

RoutingStats &
RoutingExperiment::GetRoutingStats ()
{
  return routingStats;
}
/*
void RoutingExperiment::ComputePir (Ptr<const Packet> packet, const Address &from)
{
  NS_UNUSED (from);
  if (pirCounter == 0 && lastPktRxTime.GetSeconds () == 0.0)
    {
      //this the first packet, just store the time and get out
      lastPktRxTime = Simulator::Now ();
      return;
    }
  pir = pir + (Simulator::Now () - lastPktRxTime);
  lastPktRxTime = Simulator::Now ();
  pirCounter++;
}
*/
void
RoutingExperiment::CheckThroughput ()
{


  uint32_t bytesTotal = GetRoutingStats ().GetRxBytes ();
  double kbs = (bytesTotal * 8.0) / 1000;
  uint32_t packetsReceived = GetRoutingStats().GetRxPkts();

   //double PDR = (bytesTotal * 100.0) / (txByteCounter);
  //double PDR =  packetsReceived /  (txPktCounter);
  double PDR =0.0;
  uint32_t PktsSent = GetRoutingStats().GetTxPkts();
  //int wavePktsReceived = m_waveBsmHelper.GetWaveBsmStats ()->GetRxPktCount ();
 if (packetsReceived > 0)
 {
     PDR = (double) packetsReceived / (double) PktsSent;
     if (PDR > 1)
     {
       PDR = 1;
     }

 }


  // bytesTotal = 0;

  std::ofstream out (m_CSVfileName.c_str (), std::ios::app);

  out << (Simulator::Now ()).GetSeconds () << ","
      << kbs << ","
      //<< pir << ","
      << PktsSent << ","
      << packetsReceived << ","
      << PDR << ","
      << m_nSinks << ","
      << m_protocolName << ","
      << m_txp << ","
      //<< m_numVehiclesPerLane << ","
      << std::endl;


  out.close ();
  //PktsSent = 0;
  //packetsReceived = 0;
  Simulator::Schedule (Seconds (1.0), &RoutingExperiment::CheckThroughput, this);
}

Ptr<Socket>
RoutingExperiment::SetupPacketReceive (Ipv4Address addr, Ptr<Node> node)
{
  TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
  Ptr<Socket> sink = Socket::CreateSocket (node, tid);
  InetSocketAddress local = InetSocketAddress (addr, port);
  // Bind to the local address
  sink->Bind (local);
  sink->SetRecvCallback (MakeCallback (&RoutingExperiment::ReceivePacket, this));
  //sink->SetRecvCallback (MakeCallback (&ComputePir, this));
  return sink;
}

std::string
RoutingExperiment::CommandSetup (int argc, char **argv)
{
  CommandLine cmd (__FILE__);
  cmd.AddValue ("CSVfileName", "The name of the CSV output file name", m_CSVfileName);
  cmd.AddValue ("numVehiclesPerLane", "The number of vehicles per lane", m_numVehiclesPerLane);
  cmd.AddValue ("traceMobility", "Enable mobility tracing", m_traceMobility);
  //cmd.AddValue ("nSinks", "number of sinks", m_nSinks);
  cmd.AddValue ("protocol", "1=OLSR;2=AODV;3=DSDV;4=DSR", m_protocol);
  cmd.Parse (argc, argv);
  return m_CSVfileName;

  //cmd.AddValue ("interUeDistance",
                //"The distance among the UEs in the topology",
                //interUeDistance);
}

uint32_t port=9;

int
main (int argc, char *argv[])
{
  RoutingExperiment experiment;
  std::string CSVfileName = experiment.CommandSetup (argc,argv);
//  uint16_t numVehiclesPerLane = 2;
  //CommandLine cmd;
  //cmd.AddValue ("numVehiclesPerLane",
            //    "Number of vehicles per lane", numVehiclesPerLane);
  //blank out the last output file and write the column headers
  std::ofstream out (CSVfileName.c_str ());
  out << "SimulationSecond," <<
  "ReceiveRate," <<
  //"pir," <<
  "PktsSent," <<
  "PacketsReceived," <<
  "PDR,"<<
  "NumberOfSinks," <<
  "RoutingProtocol," <<
  "TransmissionPower" <<
  std::endl;
  out.close ();

  int nSinks = 2;
  double txPower = 23;
  //uint16_t interUeDistance = 44;

  //double txPower = 10;
  //double txp = 7.5;

  experiment.Run (nSinks, txPower, CSVfileName);
}

  void
  RoutingExperiment::Run (int nSinks, double txPower, std::string CSVfileName)
  {
    Packet::EnablePrinting ();
    m_nSinks = nSinks;
    m_txp = txPower;
    m_CSVfileName = CSVfileName;
  //  m_numVehiclesPerLane = numVehiclesPerLane;
    //pdr = pDr;


   //int ueNum = 50;

    //double TotalTime = 200.0;
  //  std::string rate ("2048bps");
    //std::string phyMode ("DsssRate11Mbps");
    std::string tr_name ("test.output");
    int nodeSpeed = 0; //in m/s
    int nodePause = 0; //in s
    m_protocolName = "protocol";

    uint32_t SentPackets = 0;
    uint32_t ReceivedPackets = 0;
    uint32_t LostPackets = 0;

    uint16_t interUeDistance =2;

    bool logging = false;
    std::string m_traceFile;
     //bool enableOneTxPerLane = true;
     uint16_t interLaneDist = 2;
     uint16_t interVehicleDist = 7  ;
     uint16_t numVehiclesPerLane = m_numVehiclesPerLane;
     uint16_t numLanes = 4;
     //double totalTime = 100;
     //bool printRoutes;
    // Traffic parameters (that we will use inside this script:)
    //bool useIPv6 = false; // default IPV4
    uint32_t udpPacketSizeBe = 300;
    double dataRateBe = 24; //16 kilobits per second

    // Simulation parameters.
    Time simTime = Seconds (10);
    //Sidelink bearers activation time
    Time slBearersActivationTime = Seconds (2.0);

    double speed = 5.55;
    // NR parameters. We will take the input from the command line, and then we
    // will pass them inside the NR module.
    uint16_t numerologyBwpSl = 0;
    double centralFrequencyBandSl = 5.89e9; // band n47  TDD //Here band is analogous to channel
    uint16_t bandwidthBandSl = 100; //Multiple of 100 KHz; 400 = 40 MHz
  //  double txPower = 23; //dBm
    double slProbResourceKeep = 0;
    // Where we will store the output files.
    std::string simTag = "default";
    std::string outputDir = "./";

    //Config::SetDefault ("ns3::OnOffApplication::PacketSize", UintegerValue (64));
    //Config::SetDefault ("ns3::OnOffApplication::DataRate", StringValue ("2048bps"));


    //std::cout << "TOTO" << std::endl;

    // Final simulation time is the addition of:
    //simTime + slBearersActivationTime + 10 ms additional delay for UEs to activate the bearers
    Time finalSlBearersActivationTime = slBearersActivationTime + Seconds (0.01);
    Time finalSimTime = simTime + finalSlBearersActivationTime;
    std::cout << "Final Simulation duration " << finalSimTime.GetSeconds () << std::endl;

    /*
     * Check if the frequency is in the allowed range.
     * If you need to add other checks, here is the best position to put them.
     */
    NS_ABORT_IF (centralFrequencyBandSl > 6e9);

    /*
     * If the logging variable is set to true, enable the log of some components
     * through the code. The same effect can be obtained through the use
     * of the NS_LOG environment variable:
     *
     * export NS_LOG="UdpClient=level_info|prefix_time|prefix_func|prefix_node:UdpServer=..."
     *
     * Usually, the environment variable way is preferred, as it is more customizable,
     * and more expressive.
     */
    if (logging)
      {
        LogLevel logLevel = (LogLevel)(LOG_PREFIX_FUNC | LOG_PREFIX_TIME | LOG_PREFIX_NODE | LOG_LEVEL_ALL);
        LogComponentEnable ("UdpClient", logLevel);
        LogComponentEnable ("UdpServer", logLevel);
        LogComponentEnable ("LtePdcp", logLevel);
        LogComponentEnable ("NrSlHelper", logLevel);
        LogComponentEnable ("NrSlUeRrc", logLevel);
        LogComponentEnable ("NrUePhy", logLevel);
        LogComponentEnable ("NrSpectrumPhy", logLevel);

      }

    /*
     * Default values for the simulation. We are progressively removing all
     * the instances of SetDefault, but we need it for legacy code (LTE)
     */
    Config::SetDefault ("ns3::LteRlcUm::MaxTxBufferSize", UintegerValue(999999999));
  /*
   * Create a NodeContainer for the UEs, name it as per their traffic type.
   */
  //NodeContainer ueVoiceContainer;
  //uint16_t ueNum = 3;
  //ueVoiceContainer.Create (ueNum);
 //Add for olsr
   //NS_LOG_INFO ("Create nodes.");
 // NodeContainer c;
  //c.Create (3);
  //NodeContainer n12 = NodeContainer (ueVoiceContainer.Get (0), ueVoiceContainer.Get (1));

  //NodeContainer n32 = NodeContainer (ueVoiceContainer.Get (2), ueVoiceContainer.Get (1));


  /*
   * Assign mobility to the UEs.
   *  1. Set mobility model type.
   *  2. Assign position to the UEss
   *  3. Install mobility model
   */
/*
  MobilityHelper mobility;
  mobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  Ptr<ListPositionAllocator> positionAllocUe = CreateObject<ListPositionAllocator> ();
  for (uint16_t i = 0; i < ueNum; i++)
    {
    positionAllocUe->Add (Vector (interUeDistance * i, 0.0, 0.0));
   }
//  positionAllocUe->Add (Vector (0.0, 0.0, 0.0));
//  positionAllocUe->Add (Vector (50.0, 0.0, 0.0));
//  positionAllocUe->Add (Vector (396.0, 0.0, 0.0));
  mobility.SetPositionAllocator (positionAllocUe);
  mobility.Install (ueVoiceContainer);
*/
  NodeContainer ueVoiceContainer;
  ueVoiceContainer.Create (numVehiclesPerLane * numLanes);

  std::cout << "Total UEs = " << ueVoiceContainer.GetN () << std::endl;

  Ptr<GridPositionAllocator>  gridPositionAllocator;
  gridPositionAllocator = CreateObject<GridPositionAllocator> ();
  gridPositionAllocator->SetAttribute ("GridWidth", UintegerValue (numVehiclesPerLane));
  gridPositionAllocator->SetAttribute ("MinX", DoubleValue (0.0));
  gridPositionAllocator->SetAttribute ("MinY", DoubleValue (0.0));
  gridPositionAllocator->SetAttribute ("Z", DoubleValue (1.6));
  gridPositionAllocator->SetAttribute ("DeltaX", DoubleValue (interVehicleDist));
  gridPositionAllocator->SetAttribute ("DeltaY", DoubleValue (interLaneDist));
  gridPositionAllocator->SetAttribute ("LayoutType", EnumValue (GridPositionAllocator::ROW_FIRST));

  MobilityHelper mobility;
  mobility.SetMobilityModel ("ns3::ConstantVelocityMobilityModel");
  mobility.SetPositionAllocator (gridPositionAllocator);

  mobility.Install (ueVoiceContainer);

  for (int i = 0; i < numLanes * numVehiclesPerLane; i++)
    {
      ueVoiceContainer.Get (i)->GetObject<ConstantVelocityMobilityModel> ()->SetVelocity (Vector (speed, 0.0, 0.0));
    }

  //return ueVoiceContainer;
  //Use a trace file with ns2MobilityHelper
//  m_traceFile = "src/wave/examples/low99-ct-unterstrass-1day.filt.7.adj.mob";
   //Ns2MobilityHelper ns2 = Ns2MobilityHelper (m_traceFile);
   //ns2.Install ();
  /*
    MobilityHelper mobility;
  int64_t streamIndex = 0; // used to get consistent mobility across scenarios

  ObjectFactory pos;
  pos.SetTypeId ("ns3::RandomRectanglePositionAllocator");
  pos.Set ("X", StringValue ("ns3::UniformRandomVariable[Min=0.0|Max=300.0]"));
  pos.Set ("Y", StringValue ("ns3::UniformRandomVariable[Min=0.0|Max=1500.0]"));

  Ptr<PositionAllocator> taPositionAlloc = pos.Create ()->GetObject<PositionAllocator> ();
  streamIndex += taPositionAlloc->AssignStreams (streamIndex);

  std::stringstream ssSpeed;
  ssSpeed << "ns3::UniformRandomVariable[Min=0.0|Max=" << nodeSpeed << "]";
  std::stringstream ssPause;
  ssPause << "ns3::ConstantRandomVariable[Constant=" << nodePause << "]";
  mobility.SetMobilityModel ("ns3::RandomWaypointMobilityModel",
                                  "Speed", StringValue (ssSpeed.str ()),
                                  "Pause", StringValue (ssPause.str ()),
                                  "PositionAllocator", PointerValue (taPositionAlloc));
  mobility.SetPositionAllocator (taPositionAlloc);
  mobility.Install (ueVoiceContainer);
  */
  //streamIndex += mobility.AssignStreams (ueVoiceContainer, streamIndex);
  //NS_UNUSED (streamIndex); // From this point, streamIndex is unused
/*
 * Setup the NR module. We create the various helpers needed for the
 * NR simulation:
 * - EpcHelper, which will setup the core network
 * - NrHelper, which takes care of creating and connecting the various
 * part of the NR stack
 */
  Ptr<NrPointToPointEpcHelper> epcHelper = CreateObject<NrPointToPointEpcHelper> ();
  Ptr<NrHelper> nrHelper = CreateObject<NrHelper> ();

  // Put the pointers inside nrHelper
  nrHelper->SetEpcHelper (epcHelper);
  //nrHelper->SetPathlossModelType("ns3::RangePropagationLossModel");

  /*
   * Spectrum division. We create one operational band, containing
   * one component carrier, and a single bandwidth part
   * centered at the frequency specified by the input parameters.
   * We will use the StreetCanyon channel modeling.
   */
  BandwidthPartInfoPtrVector allBwps;
  CcBwpCreator ccBwpCreator;
  const uint8_t numCcPerBand = 1;

  /* Create the configuration for the CcBwpHelper. SimpleOperationBandConf
   * creates a single BWP per CC
   */
  CcBwpCreator::SimpleOperationBandConf bandConfSl (centralFrequencyBandSl, bandwidthBandSl, numCcPerBand, BandwidthPartInfo::V2V_Highway);

  //CcBwpCreator::SimpleOperationBandConf bandConfSl (centralFrequencyBandSl, bandwidthBandSl, numCcPerBand, BandwidthPartInfo::V2V_Urban);

  // By using the configuration created, it is time to make the operation bands
  OperationBandInfo bandSl = ccBwpCreator.CreateOperationBandContiguousCc (bandConfSl);

  /*
   * The configured spectrum division is:
   * ------------Band1--------------
   * ------------CC1----------------
   * ------------BwpSl--------------
   */

  /*
   * Attributes of ThreeGppChannelModel still cannot be set in our way.
   * TODO: Coordinate with Tommaso
   */
  Config::SetDefault ("ns3::ThreeGppChannelModel::UpdatePeriod",TimeValue (MilliSeconds(100)));
  nrHelper->SetChannelConditionModelAttribute ("UpdatePeriod", TimeValue (MilliSeconds (0)));
  nrHelper->SetPathlossAttribute ("ShadowingEnabled", BooleanValue (false));
  //Config::SetDefault ("ns3::RangePropagationLossModel::MaxRange", BooleanValue (150));
  //nrHelper->SetPathlossAttribute ("MaxRange", BooleanValue (100));
  //nrHelper->SetPathlossAttribute ("MaxRange", BooleanValue (150));
  //nrHelper->SetPathlossModelType("ns3::RangePropagationLossModel");

  /*
   * Initialize channel and pathloss, plus other things inside bandSl. If needed,
   * the band configuration can be done manually, but we leave it for more
   * sophisticated examples. For the moment, this method will take care
   * of all the spectrum initialization needs.
   */
  nrHelper->InitializeOperationBand (&bandSl);
  allBwps = CcBwpCreator::GetAllBwps ({bandSl});
  /*
   * allBwps contains all the spectrum configuration needed for the nrHelper.
   *
   * Now, we can setup the attributes. We can have three kind of attributes:
   * (i) parameters that are valid for all the bandwidth parts and applies to
   * all nodes, (ii) parameters that are valid for all the bandwidth parts
   * and applies to some node only, and (iii) parameters that are different for
   * every bandwidth parts. The approach is:
   *
   * - for (i): Configure the attribute through the helper, and then install;
   * - for (ii): Configure the attribute through the helper, and then install
   * for the first set of nodes. Then, change the attribute through the helper,
   * and install again;
   * - for (iii): Install, and then configure the attributes by retrieving
   * the pointer needed, and calling "SetAttribute" on top of such pointer.
   *
   */

  Packet::EnableChecking ();
  Packet::EnablePrinting ();

  /*
   *  Case (i): Attributes valid for all the nodes
   */
  // Core latency
  epcHelper->SetAttribute ("S1uLinkDelay", TimeValue (MilliSeconds (0)));

  /*
   * Antennas for all the UEs
   * We are not using beamforming in SL, rather we are using
   * quasi-omnidirectional transmission and reception, which is the default
   * configuration of the beams.
   */
  nrHelper->SetUeAntennaAttribute ("NumRows", UintegerValue (1));
  nrHelper->SetUeAntennaAttribute ("NumColumns", UintegerValue (2));
  nrHelper->SetUeAntennaAttribute ("AntennaElement", PointerValue (CreateObject<IsotropicAntennaModel> ()));

  nrHelper->SetUePhyAttribute ("TxPower", DoubleValue (txPower));

  //NR Sidelink attribute of UE MAC, which are would be common for all the UEs
  nrHelper->SetUeMacAttribute ("EnableSensing", BooleanValue (true));
  nrHelper->SetUeMacAttribute ("T1", UintegerValue (2));
  nrHelper->SetUeMacAttribute ("T2", UintegerValue (33));
  nrHelper->SetUeMacAttribute ("ActivePoolId", UintegerValue (0));
  nrHelper->SetUeMacAttribute ("ReservationPeriod", TimeValue (MilliSeconds(100)));
  nrHelper->SetUeMacAttribute ("NumSidelinkProcess", UintegerValue (4));
  nrHelper->SetUeMacAttribute ("EnableBlindReTx", BooleanValue (true));


  uint8_t bwpIdForGbrMcptt = 0;

  nrHelper->SetBwpManagerTypeId (TypeId::LookupByName ("ns3::NrSlBwpManagerUe"));
  //following parameter has no impact at the moment because:
  //1. No support for PQI based mapping between the application and the LCs
  //2. No scheduler to consider PQI
  //However, till such time all the NR SL examples should use GBR_MC_PUSH_TO_TALK
  //because we hard coded the PQI 65 in UE RRC.
  nrHelper->SetUeBwpManagerAlgorithmAttribute ("GBR_MC_PUSH_TO_TALK", UintegerValue (bwpIdForGbrMcptt));

  std::set<uint8_t> bwpIdContainer;
  bwpIdContainer.insert (bwpIdForGbrMcptt);

  /*
   * We have configured the attributes we needed. Now, install and get the pointers
   * to the NetDevices, which contains all the NR stack:
   */
  NetDeviceContainer ueVoiceNetDev = nrHelper->InstallUeDevice (ueVoiceContainer, allBwps);

  // When all the configuration is done, explicitly call UpdateConfig ()
  for (auto it = ueVoiceNetDev.Begin (); it != ueVoiceNetDev.End (); ++it)
    {
      DynamicCast<NrUeNetDevice> (*it)->UpdateConfig ();
    }

    /*
     * Configure Sidelink. We create the following helpers needed for the
     * NR Sidelink, i.e., V2X simulation:
     * - NrSlHelper, which will configure the UEs protocol stack to be ready to
     *   perform Sidelink related procedures.
     * - EpcHelper, which takes care of triggering the call to EpcUeNas class
     *   to establish the NR Sidelink bearer (s). We note that, at this stage
     *   just communicate the pointer of already instantiated EpcHelper object,
     *   which is the same pointer communicated to the NrHelper above.
     */
    Ptr<NrSlHelper> nrSlHelper = CreateObject <NrSlHelper> ();
    // Put the pointers inside NrSlHelper
    nrSlHelper->SetEpcHelper (epcHelper);

    /*
     * Set the SL error model and AMC
     * Error model type: ns3::NrEesmCcT1, ns3::NrEesmCcT2, ns3::NrEesmIrT1,
     *                   ns3::NrEesmIrT2, ns3::NrLteMiErrorModel
     * AMC type: NrAmc::ShannonModel or NrAmc::ErrorModel
     */
    std::string errorModel = "ns3::NrEesmIrT1";
    nrSlHelper->SetSlErrorModel (errorModel);
    nrSlHelper->SetUeSlAmcAttribute ("AmcModel", EnumValue (NrAmc::ErrorModel));

    /*
     * Set the SL scheduler attributes
     * In this example we use NrSlUeMacSchedulerSimple scheduler, which uses
     * fix MCS value
     */
    nrSlHelper->SetNrSlSchedulerTypeId (NrSlUeMacSchedulerSimple::GetTypeId());
    nrSlHelper->SetUeSlSchedulerAttribute ("FixNrSlMcs", BooleanValue (true));
    nrSlHelper->SetUeSlSchedulerAttribute ("InitialNrSlMcs", UintegerValue (14));
    /*
     * Very important method to configure UE protocol stack, i.e., it would
     * configure all the SAPs among the layers, setup callbacks, configure
     * error model, configure AMC, and configure ChunkProcessor in Interference
     * API.
     */
    nrSlHelper->PrepareUeForSidelink (ueVoiceNetDev, bwpIdContainer);

    /*
     * Start preparing for all the sub Structs/RRC Information Element (IEs)
     * of LteRrcSap::SidelinkPreconfigNr. This is the main structure, which would
     * hold all the pre-configuration related to Sidelink.
     */

    //SlResourcePoolNr IE
    LteRrcSap::SlResourcePoolNr slResourcePoolNr;
    //get it from pool factory
    Ptr<NrSlCommPreconfigResourcePoolFactory> ptrFactory = Create<NrSlCommPreconfigResourcePoolFactory> ();

    /*
     * Above pool factory is created to help the users of the simulator to create
     * a pool with valid default configuration. Please have a look at the
     * constructor of NrSlCommPreconfigResourcePoolFactory class.
     *
     * In the following, we show how one could change those default pool parameter
     * values as per the need.
     */
    std::vector <std::bitset<1> > slBitmap = {1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1};
    ptrFactory->SetSlTimeResources (slBitmap);
    ptrFactory->SetSlSensingWindow (1100); // T0 in ms
    ptrFactory->SetSlSelectionWindow (5);
    ptrFactory->SetSlFreqResourcePscch (10); // PSCCH RBs
    ptrFactory->SetSlSubchannelSize (50);
    ptrFactory->SetSlMaxNumPerReserve (3);
    //Once parameters are configured, we can create the pool
    LteRrcSap::SlResourcePoolNr pool = ptrFactory->CreatePool ();
    slResourcePoolNr = pool;

    //Configure the SlResourcePoolConfigNr IE, which hold a pool and its id
    LteRrcSap::SlResourcePoolConfigNr slresoPoolConfigNr;
    slresoPoolConfigNr.haveSlResourcePoolConfigNr = true;
    //Pool id, ranges from 0 to 15
    uint16_t poolId = 0;
    LteRrcSap::SlResourcePoolIdNr slResourcePoolIdNr;
    slResourcePoolIdNr.id = poolId;
    slresoPoolConfigNr.slResourcePoolId = slResourcePoolIdNr;
    slresoPoolConfigNr.slResourcePool = slResourcePoolNr;

    //Configure the SlBwpPoolConfigCommonNr IE, which hold an array of pools
    LteRrcSap::SlBwpPoolConfigCommonNr slBwpPoolConfigCommonNr;
    //Array for pools, we insert the pool in the array as per its poolId
    slBwpPoolConfigCommonNr.slTxPoolSelectedNormal [slResourcePoolIdNr.id] = slresoPoolConfigNr;

    //Configure the BWP IE
    LteRrcSap::Bwp bwp;
    bwp.numerology = numerologyBwpSl;
    bwp.symbolsPerSlots = 14;
    bwp.rbPerRbg = 1;
    bwp.bandwidth = bandwidthBandSl;

    //Configure the SlBwpGeneric IE
    LteRrcSap::SlBwpGeneric slBwpGeneric;
    slBwpGeneric.bwp = bwp;
    slBwpGeneric.slLengthSymbols = LteRrcSap::GetSlLengthSymbolsEnum (14);
    slBwpGeneric.slStartSymbol = LteRrcSap::GetSlStartSymbolEnum (0);

    //Configure the SlBwpConfigCommonNr IE
    LteRrcSap::SlBwpConfigCommonNr slBwpConfigCommonNr;
    slBwpConfigCommonNr.haveSlBwpGeneric = true;
    slBwpConfigCommonNr.slBwpGeneric = slBwpGeneric;
    slBwpConfigCommonNr.haveSlBwpPoolConfigCommonNr = true;
    slBwpConfigCommonNr.slBwpPoolConfigCommonNr = slBwpPoolConfigCommonNr;

    //Configure the SlFreqConfigCommonNr IE, which hold the array to store
    //the configuration of all Sidelink BWP (s).
    LteRrcSap::SlFreqConfigCommonNr slFreConfigCommonNr;
    //Array for BWPs. Here we will iterate over the BWPs, which
    //we want to use for SL.
    for (const auto &it:bwpIdContainer)
      {
        // it is the BWP id
        slFreConfigCommonNr.slBwpList [it] = slBwpConfigCommonNr;
      }

      //Configure the TddUlDlConfigCommon IE
      LteRrcSap::TddUlDlConfigCommon tddUlDlConfigCommon;
      tddUlDlConfigCommon.tddPattern = "DL|DL|DL|F|UL|UL|UL|UL|UL|UL|";

      //Configure the SlPreconfigGeneralNr IE
      LteRrcSap::SlPreconfigGeneralNr slPreconfigGeneralNr;
      slPreconfigGeneralNr.slTddConfig = tddUlDlConfigCommon;

      //Configure the SlUeSelectedConfig IE
      LteRrcSap::SlUeSelectedConfig slUeSelectedPreConfig;
      slUeSelectedPreConfig.slProbResourceKeep = 0;
      //Configure the SlPsschTxParameters IE
      LteRrcSap::SlPsschTxParameters psschParams;
      psschParams.slMaxTxTransNumPssch = 5;
      //Configure the SlPsschTxConfigList IE
      LteRrcSap::SlPsschTxConfigList pscchTxConfigList;
      pscchTxConfigList.slPsschTxParameters [0] = psschParams;
      slUeSelectedPreConfig.slPsschTxConfigList = pscchTxConfigList;

      /*
       * Finally, configure the SidelinkPreconfigNr This is the main structure
       * that needs to be communicated to NrSlUeRrc class
       */
      LteRrcSap::SidelinkPreconfigNr slPreConfigNr;
      slPreConfigNr.slPreconfigGeneral = slPreconfigGeneralNr;
      slPreConfigNr.slUeSelectedPreConfig = slUeSelectedPreConfig;
      slPreConfigNr.slPreconfigFreqInfoList [0] = slFreConfigCommonNr;

      //Communicate the above pre-configuration to the NrSlHelper
      nrSlHelper->InstallNrSlPreConfiguration (ueVoiceNetDev, slPreConfigNr);




      /****************************** End SL Configuration ***********************/

      /*
       * Fix the random streams
       */
      int64_t stream = 1;
      stream += nrHelper->AssignStreams (ueVoiceNetDev, stream);
      stream += nrSlHelper->AssignStreams (ueVoiceNetDev, stream);

      AodvHelper aodv;
      OlsrHelper olsr;
      Ipv4StaticRoutingHelper staticRouting;

      Ipv4ListRoutingHelper list;
      list.Add (staticRouting, 0);
      //list.Add (olsr, 10);

      InternetStackHelper internet;
      //internet.SetRoutingHelper (list); // has effect on the next Install ()
    //  internet.Install (ueVoiceContainer);


        switch (m_protocol)
          {
         case 1:
         list.Add (olsr, 100);
         m_protocolName = "OLSR";
         break;
          case 2:
           list.Add (aodv, 100);
           m_protocolName = "AODV";
           break;
          //case 3:
            //list.Add (dsdv, 100);
            //m_protocolName = "DSDV";
          //  break;
          //case 4:
            //m_protocolName = "DSR";
            //break;
          default:
            NS_FATAL_ERROR ("No such protocol:" << m_protocol);
          }

        if (m_protocol < 4)
          {
            internet.SetRoutingHelper (list);
            internet.Install (ueVoiceContainer);
          }

       /*
          NodeContainer txSlUes;
          NodeContainer rxSlUes;
          NetDeviceContainer txSlUesNetDevice;
          NetDeviceContainer rxSlUesNetDevice;

          //if (enableOneTxPerLane)
           // {
              for (uint16_t i = 1; i < ueVoiceContainer.GetN (); i++)
                {
                  //for each lane one transmitter
                  //if (i % numVehiclesPerLane == 0)
                    //{
                     // uint16_t firstIndex = (i - numVehiclesPerLane) + 1;
                      //uint16_t txNodeId = (firstIndex + i) / 2;
                     // txNodeId = txNodeId - 1; //node id starts from 0
                      txSlUes.Add (ueVoiceContainer.Get (i));
                      Ptr<NetDevice> dev = ueVoiceContainer.Get (i)->GetDevice (0);
                      Ptr<NrUeNetDevice> ueDev = DynamicCast <NrUeNetDevice> (dev);
                      NS_ABORT_MSG_IF (ueDev == nullptr, "Device 0 is not the NrUeNetDevice");

                      txSlUesNetDevice.Add (ueVoiceContainer.Get (i)->GetDevice (0));

                   // }
                }

              //all node ids, which are not in txSlUes container are Rx node ids
              for (uint16_t i = 0; i  < ueVoiceContainer.GetN (); i++)
                {
                  Ptr<Node> node = ueVoiceContainer.Get (i);
                  if (!txSlUes.Contains (node->GetId ()))
                    {
                      rxSlUes.Add (node);
                      rxSlUesNetDevice.Add (node->GetDevice (0));
                      Ptr<NetDevice> dev = ueVoiceContainer.Get (node->GetId ())->GetDevice (0);
                      Ptr<NrUeNetDevice> ueDev = DynamicCast <NrUeNetDevice> (dev);
                      NS_ABORT_MSG_IF (ueDev == nullptr, "Device 0 is not the NrUeNetDevice");
                    }
                  else
                   {
                     txSlUes.Add (ueVoiceContainer);
                     rxSlUes.Add (ueVoiceContainer);
                     txSlUesNetDevice.Add (ueVoiceNetDev);
                     rxSlUesNetDevice.Add (ueVoiceNetDev);
                    }
                  }
*/
   //Ipv4AddressHelper groupAddress4;
   //groupAddress4.SetBase ("10.1.1.0", "255.255.255.0");

   Ipv4Address groupAddress4 ("255.255.255.0");
   Address remoteAddress;
   Address localAddress;
   remoteAddress = InetSocketAddress (groupAddress4, port);
   localAddress = InetSocketAddress (Ipv4Address::GetAny (), port);
   Ipv4InterfaceContainer ueIpIface;
   ueIpIface = epcHelper->AssignUeIpv4Address (ueVoiceNetDev);
   uint32_t dstL2Id = 255;
  Ptr<LteSlTft> tft;

          tft = Create<LteSlTft> (LteSlTft::Direction::BIDIRECTIONAL, LteSlTft::CommType::Broadcast, groupAddress4, dstL2Id);
          //Set Sidelink bearers
          nrSlHelper->ActivateNrSlBearer (finalSlBearersActivationTime, ueVoiceNetDev, tft);

   /*
   Ipv4InterfaceContainer ueIpIface;
   ueIpIface = epcHelper->AssignUeIpv4Address (ueVoiceNetDev);
          uint32_t dstL2Id = 255;
          Ipv4Address groupAddress4 ("225.0.0.0");     //use multicast address as destination
          Ipv6Address groupAddress6 ("ff0e::1");     //use multicast address as destination
          Address remoteAddress;
          Address localAddress;
          uint16_t port = 8000;
          Ptr<LteSlTft> tft;
          if (!useIPv6)
            {
              Ipv4InterfaceContainer ueIpIface;
              ueIpIface = epcHelper->AssignUeIpv4Address (ueVoiceNetDev);

              // set the default gateway for the UE
              Ipv4StaticRoutingHelper ipv4RoutingHelper;
              for (uint32_t u = 0; u < ueVoiceContainer.GetN (); ++u)
                {
                  Ptr<Node> ueNode = ueVoiceContainer.Get (u);
                  // Set the default gateway for the UE
                  Ptr<Ipv4StaticRouting> ueStaticRouting = ipv4RoutingHelper.GetStaticRouting (ueNode->GetObject<Ipv4> ());
                  ueStaticRouting->SetDefaultRoute (epcHelper->GetUeDefaultGatewayAddress (), 1);
                }
              remoteAddress = InetSocketAddress (groupAddress4, port);
              localAddress = InetSocketAddress (Ipv4Address::GetAny (), port);

              tft = Create<LteSlTft> (LteSlTft::Direction::TRANSMIT, LteSlTft::CommType::GroupCast, groupAddress4, dstL2Id);
              //Set Sidelink bearers
              nrSlHelper->ActivateNrSlBearer (slBearersActivationTime, txSlUesNetDevice, tft);

              tft = Create<LteSlTft> (LteSlTft::Direction::RECEIVE, LteSlTft::CommType::GroupCast, groupAddress4, dstL2Id);
              //Set Sidelink bearers
              nrSlHelper->ActivateNrSlBearer (slBearersActivationTime, rxSlUesNetDevice, tft);
            }
          else
            {
              Ipv6InterfaceContainer ueIpIface;
              ueIpIface = epcHelper->AssignUeIpv6Address (ueVoiceNetDev);

              // set the default gateway for the UE
              Ipv6StaticRoutingHelper ipv6RoutingHelper;
              for (uint32_t u = 0; u < ueVoiceContainer.GetN (); ++u)
                {
                  Ptr<Node> ueNode = ueVoiceContainer.Get (u);
                  // Set the default gateway for the UE
                  Ptr<Ipv6StaticRouting> ueStaticRouting = ipv6RoutingHelper.GetStaticRouting (ueNode->GetObject<Ipv6> ());
                  ueStaticRouting->SetDefaultRoute (epcHelper->GetUeDefaultGatewayAddress6 (), 1);
                }
              remoteAddress = Inet6SocketAddress (groupAddress6, port);
              localAddress = Inet6SocketAddress (Ipv6Address::GetAny (), port);

              tft = Create<LteSlTft> (LteSlTft::Direction::TRANSMIT, LteSlTft::CommType::GroupCast, groupAddress4, dstL2Id);
              //Set Sidelink bearers for transmitting UEs
              nrSlHelper->ActivateNrSlBearer (slBearersActivationTime, txSlUesNetDevice, tft);

              tft = Create<LteSlTft> (LteSlTft::Direction::RECEIVE, LteSlTft::CommType::GroupCast, groupAddress4, dstL2Id);
              //Set Sidelink bearers for receiving UEs
              nrSlHelper->ActivateNrSlBearer (slBearersActivationTime, rxSlUesNetDevice, tft);
            }
      */

    Ptr<UniformRandomVariable> startTimeSeconds = CreateObject<UniformRandomVariable> ();
    startTimeSeconds->SetStream (stream);
    startTimeSeconds->SetAttribute ("Min", DoubleValue (0));
    startTimeSeconds->SetAttribute ("Max", DoubleValue (12));

   OnOffHelper sidelinkClient ("ns3::UdpSocketFactory", Address ());
   //sidelinkClient.SetAttribute ("EnableSeqTsSizeHeader", BooleanValue (true));
  sidelinkClient.SetAttribute ("OnTime", StringValue ("ns3::ConstantRandomVariable[Constant=1.0]"));
  sidelinkClient.SetAttribute ("OffTime", StringValue ("ns3::ConstantRandomVariable[Constant=0.0]"));

  //std::string dataRateBeString  = std::to_string (dataRateBe) + "kb/s";
  //std::cout << "Data rate " << DataRate (dataRateBeString) << std::endl;
  //sidelinkClient.SetConstantRate (DataRate (dataRateBe), udpPacketSizeBe);
   double realAppStart = 0.0;
   double txAppDuration = 0.0;
   //Output app start, stop and duration
   //double realAppStart =  finalSlBearersActivationTime.GetSeconds () + ((double)udpPacketSizeBe * 8.0 / (DataRate (dataRateBeString).GetBitRate ()));
 //realAppStart =  finalSlBearersActivationTime.GetSeconds () + ((double)udpPacketSizeBe * 8.0 / (DataRate (dataRateBe).GetBitRate ()));
 realAppStart = finalSlBearersActivationTime.GetSeconds () ;
 double appStopTime = (finalSimTime).GetSeconds ();
  txAppDuration = appStopTime - realAppStart;

   std::cout << "App start time " << realAppStart << " sec" << std::endl;
   std::cout << "App stop time " << appStopTime << " sec" << std::endl;
   std::cout << "Tx App duration " << std::defaultfloat << txAppDuration  << " sec"  << std::endl;



          //ApplicationContainer clientApps;
        //  double realAppStart = 0.0;
          //double realAppStopTime = 0.0;
          //double txAppDuration = 0.0;

       /*
          for (uint16_t i = 0; i  < nSinks; i++)
            {
              if (m_protocol != 0)
                {

                  Ptr<Socket> sink = SetupPacketReceive (ueIpIface.GetAddress (i), ueVoiceContainer.Get (i));
                }

              AddressValue remoteAddress (InetSocketAddress (ueIpIface.GetAddress (i), port));
              sidelinkClient.SetAttribute ("Remote", remoteAddress);
              clientApps.Add (sidelinkClient.Install (txSlUes.Get (i + nSinks)));
              //double jitter = startTimeSeconds->GetValue ();
              Time appStart = slBearersActivationTime;
              clientApps.Get (i)->SetStartTime (appStart);
              //onoff application will send the first packet at :
              //slBearersActivationTime + random jitter + ((Pkt size in bits) / (Data rate in bits per sec))
              realAppStart =  slBearersActivationTime.GetSeconds () + ((double)udpPacketSizeBe * 8.0 / (DataRate (dataRateBeString).GetBitRate ()));
              //realAppStopTime = realAppStart + simTime.GetSeconds ();
              clientApps.Get (i)->SetStopTime (Seconds (appStopTime));
              //txAppDuration = realAppStopTime - realAppStart;
         }
     */
  {
    for (int  i = 0; i < nSinks; i++)
      {
        // protocol == 0 means no routing data, WAVE BSM only
        // so do not set up sink
        if (m_protocol != 0)
          {

            Ptr<Socket> sink = SetupPacketReceive (ueIpIface.GetAddress (i), ueVoiceContainer.Get (i));
          }

        AddressValue remoteAddress (InetSocketAddress (ueIpIface.GetAddress (i), port));
      //  UdpEchoClientHelper sidelinkClient (ueIpIface.GetAddress (i), port);
         //Time interval = Seconds (1.);
        //   uint32_t maxPacketCount = 1;
        //sidelinkClient.SetAttribute ("MaxPackets", UintegerValue (maxPacketCount));
        // sidelinkClient.SetAttribute ("PacketSize", UintegerValue (udpPacketSizeBe));
         //sidelinkClient.SetAttribute ("Interval", TimeValue (interval));

        sidelinkClient.SetAttribute ("Remote", remoteAddress);
        //for (uint16_t i = 0; i  < txSlUes.GetN (); i++)
        //  {
        //ApplicationContainer temp = sidelinkClient.Install (txSlUes.Get (i));
        ApplicationContainer temp = sidelinkClient.Install (ueVoiceContainer.Get (i+ nSinks));
        temp.Start (Seconds (2.0));
      }
/*std::string dataRateBeString  = std::to_string (dataRateBe) + "kb/s";
      std::cout << "Data rate " << DataRate (dataRateBeString) << std::endl;
    //  sidelinkClient.SetConstantRate (DataRate (dataRateBeString), udpPacketSizeBe);
       double realAppStart = 0.0;
       double txAppDuration = 0.0;
       //Output app start, stop and duration
       //double realAppStart =  finalSlBearersActivationTime.GetSeconds () + ((double)udpPacketSizeBe * 8.0 / (DataRate (dataRateBeString).GetBitRate ()));
     realAppStart =  finalSlBearersActivationTime.GetSeconds () + ((double)udpPacketSizeBe * 8.0 / (DataRate (dataRateBeString).GetBitRate ()));
    //realAppStart = finalSlBearersActivationTime.GetSeconds () ;
    double appStopTime = (finalSimTime).GetSeconds ();
       txAppDuration = appStopTime - realAppStart;

       std::cout << "App start time " << realAppStart << " sec" << std::endl;
       std::cout << "App stop time " << appStopTime << " sec" << std::endl;
       std::cout << "Tx App duration " << std::defaultfloat << txAppDuration  << " sec"  << std::endl;

*/
        //temp.Stop (Seconds (m_TotalSimTime));
      //}

    // Create packet sinks to receive these packets
  /*
    Ptr<Socket> sink = SetupPacketReceive (ueIpIface.GetAddress (0), ueVoiceContainer.Get (0));

    AddressValue remoteAddress (InetSocketAddress (ueIpIface.GetAddress (0), port));
    sidelinkClient.SetAttribute ("Remote", remoteAddress);
    ApplicationContainer temp = sidelinkClient.Install (ueVoiceContainer.Get (7));
    temp.Start (Seconds (2.0));
   */
/*
    ApplicationContainer serverApps;
    PacketSinkHelper sidelinkSink ("ns3::UdpSocketFactory", ueIpIface.GetAddress (0));
    sidelinkSink.SetAttribute ("EnableSeqTsSizeHeader", BooleanValue (true));
    serverApps = sidelinkSink.Install (ueVoiceContainer.Get (0));
    serverApps.Start (Seconds (2.0));
    */

  }
//std::cout << "TOTO" << std::endl;
/*  std::stringstream ss;
  ss << ueNum;
  std::string nodes = ss.str ();

  std::stringstream ss2;
  ss2 << nodeSpeed;
  std::string sNodeSpeed = ss2.str ();

  std::stringstream ss3;
  ss3 << nodePause;
  std::string sNodePause = ss3.str ();

  std::stringstream ss4;
  ss4 << dataRateBe;
  std::string sRate = ss4.str ();
*/
  /*
   * Hook the traces, to be used to compute average PIR and to data to be
   * stored in a database
   */

  //Trace receptions; use the following to be robust to node ID changes
   std::ostringstream path;
   //path << "/NodeList/" << ueVoiceContainer.Get (0)->GetId () <<"/ApplicationList/0/ns3::olsr::RoutingProtocol/Rx";
   //Config::Connect(path.str (), MakeCallback (&RoutingExperiment::ComputePir, this));
  // path.str ("");


  //std::vector<ns3::olsr::RoutingTableEntry> GetRoutingTableEntries(void) const;
  //std::cout << "TOTO" << std::endl;

  // extract the routing table of node
  /*
  Ptr<Node> node = ueVoiceContainer.Get (0);
  std::vector<ns3::olsr::RoutingTableEntry> entries = std::vector<ns3::olsr::RoutingTableEntry>();
  Ptr<ns3::olsr::RoutingProtocol> olsrOb = node->GetObject<ns3::olsr::RoutingProtocol>();
  entries = olsrOb->GetRoutingTableEntries();
  */
  /*path << "/NodeList/" << ueVoiceContainer.Get (0)->GetId () << "/$ns3::olsr::RoutingProtocol";
  Config::ConnectWithoutContext(path.str (), MakeCallback (&ns3::olsr::RoutingProtocol::TableChangeTracedCallback);
  path.str ("");*/

  //  path << "/NodeList/" << ueVoiceContainer.Get (17)->GetId () << "/ApplicationList/0/$ns3::OnOffApplication/Tx";
  //  Config::ConnectWithoutContext(path.str (), MakeCallback (&TransmitPacket));
  //  path.str ("");
  path << "/NodeList/*/ApplicationList/*/$ns3::OnOffApplication/Tx";
  Config::Connect (path.str (), MakeCallback (&RoutingExperiment::OnOffTrace, this));
  path.str ("");

  //Datebase setup
  std::string exampleName = simTag + "-" + "test-olsr";
  SQLiteOutput db (outputDir + exampleName + ".db", exampleName);

  UeMacPscchTxOutputStats pscchStats;
  pscchStats.SetDb (&db, "pscchTxUeMac");
  Config::ConnectWithoutContext ("/NodeList/*/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUeMac/SlPscchScheduling",
                                   MakeBoundCallback (&NotifySlPscchScheduling, &pscchStats));

  UeMacPsschTxOutputStats psschStats;
  psschStats.SetDb (&db, "psschTxUeMac");
  Config::ConnectWithoutContext ("/NodeList/*/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUeMac/SlPsschScheduling",
                                   MakeBoundCallback (&NotifySlPsschScheduling, &psschStats));


  UePhyPscchRxOutputStats pscchPhyStats;
  pscchPhyStats.SetDb (&db, "pscchRxUePhy");
  Config::ConnectWithoutContext ("/NodeList/*/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUePhy/SpectrumPhy/RxPscchTraceUe",
                                   MakeBoundCallback (&NotifySlPscchRx, &pscchPhyStats));

  UePhyPsschRxOutputStats psschPhyStats;
  psschPhyStats.SetDb (&db, "psschRxUePhy");
  Config::ConnectWithoutContext ("/NodeList/*/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUePhy/SpectrumPhy/RxPsschTraceUe",
                                   MakeBoundCallback (&NotifySlPsschRx, &psschPhyStats));

   UeRlcRxOutputStats ueRlcRxStats;
   ueRlcRxStats.SetDb (&db, "rlcRx");
   Config::ConnectWithoutContext ("/NodeList/*/DeviceList/*/$ns3::NrUeNetDevice/ComponentCarrierMapUe/*/NrUeMac/RxRlcPduWithTxRnti",
                                    MakeBoundCallback (&NotifySlRlcPduRx, &ueRlcRxStats));

  UeToUePktTxRxOutputStats pktStats;
  pktStats.SetDb (&db, "pktTxRx");


  // Set Tx traces
  /*
       for (uint16_t ac = 0; ac < clientApps.GetN (); ac++)
         {
           Ipv4Address localAddrs =  clientApps.Get (ac)->GetNode ()->GetObject<Ipv4L3Protocol> ()->GetAddress (1,0).GetLocal ();
           std::cout << "UE2 address: " << localAddrs << std::endl;
        //  clientApps.Get (ac)->TraceConnect ("TxWithSeqTsSize", "tx", MakeBoundCallback (&UePacketTraceDb, &pktStats, ueVoiceContainer.Get (0), localAddrs));
          clientApps.Get (ac)->TraceConnect ("TxWithSeqTsSize", "tx", MakeBoundCallback (&UePacketTraceDb, &pktStats,  clientApps.Get (ac)->GetNode (), localAddrs));

         }

       // Set Rx traces
       for (uint16_t ac = 0; ac < serverApps.GetN (); ac++)
         {
           Ipv4Address localAddrs =  serverApps.Get (ac)->GetNode ()->GetObject<Ipv4L3Protocol> ()->GetAddress (1,0).GetLocal ();
           std::cout << "UE1 address: " << localAddrs << std::endl;
           serverApps.Get (ac)->TraceConnect ("RxWithSeqTsSize", "rx", MakeBoundCallback (&UePacketTraceDb, &pktStats, ueVoiceContainer.Get (0), localAddrs));
          // serverApps.Get (ac)->TraceConnect ("RxWithSeqTsSize", "rx", MakeBoundCallback (&UePacketTraceDb, &pktStats, serverApps.Get (ac)->GetNode (), localAddrs));

         }


         //std::string tr_name ("test");
         */


         //AsciiTraceHelper ascii;
        // MobilityHelper::EnableAsciiAll (ascii.CreateFileStream (tr_name + ".mob"));
         //NrUePhy.EnablePcapAll ("test-routing");
         //PointToPoint.EnablePcapAll ("test-routing-output");
           //nrSlHelper.EnablePcapAll ("test-routing-output");
         //PointToPointHelper p2ph;
         //p2ph.EnablePcapAll("scratch/my-pcap_files");


         FlowMonitorHelper flowmon;
         Ptr<FlowMonitor> monitor = flowmon.InstallAll();

      //Time simStopTime = simTime + slBearersActivationTime + Seconds (realAppStart);

         //std::cout << "TOTO" << std::endl;
        //nrHelper->EnableTraces();
         //AsciiTraceHelper ascii;
           //nrSlHelper.EnableAsciiAll (ascii.CreateFileStream ("est-routing-output.tr"));
/*
  V2xKpi v2xKpi;
  v2xKpi.SetDbPath (outputDir + exampleName);
  v2xKpi.SetTxAppDuration (txAppDuration);
  SavePositionPerIP (&v2xKpi);
  v2xKpi.SetRangeForV2xKpis (200);

  if (generateInitialPosGnuScript)
    {
      std::string initPosFileName = "init-pos-ues-" + exampleName + ".txt";
      PrintUeInitPosToFile (initPosFileName);
    }
*/
  // Final simulation stop time is the addition of: simTime + slBearersActivationTime + realAppStart
  // realAppStart is of the last UE for which we installed the application
  //Time simStopTime = simTime + slBearersActivationTime + Seconds (realAppStart);
/*
  if (generateGifGnuScript)
    {
      std::string mobilityFileName = "mobility-" + exampleName + ".txt";
      RecordMobility (true, mobilityFileName);
      WriteGifGnuScript (mobilityFileName, finalSimTime, speed, rxSlUes.Get (0), rxSlUes.Get (rxSlUes.GetN () - 1));
    }
*/
          NS_LOG_INFO ("Run Simulation.");

          CheckThroughput ();

          Simulator::Stop (finalSimTime);
          Simulator::Run ();



         //std::cout << "Total Tx bits = " << txByteCounter * 8 << std:: endl;
         //std::cout << "Total Tx packets = " << txPktCounter << std:: endl;



         //std::cout << "Total Rx bits = " << bytesTotal* 8 << std:: endl;
         //std::cout << "Total Rx packets = " << packetsReceived << std:: endl;

        //std::cout << "Avrg thput = " << (rxByteCounter * 8) / (finalSimTime - Seconds(realAppStart)).GetSeconds () / 1000.0 << " kbps" << std:: endl;
        //std::cout << "Avrg thput = " << (bytesTotal * 8) / (finalSimTime - Seconds(realAppStart)).GetSeconds () / 1000.0 << " kbps" << std:: endl;
/*
std::cout << "TOTO" << std::endl;
        int j=0;
float AvgThroughput = 0;
Time Jitter;
Time Delay;

Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier> (flowmon.GetClassifier ());
 std::map<FlowId, FlowMonitor::FlowStats> stats = monitor->GetFlowStats ();

 for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator iter = stats.begin (); iter != stats.end (); ++iter)
   {
   Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow (iter->first);

NS_LOG_UNCOND("----Flow ID:" <<iter->first);
NS_LOG_UNCOND("Src Addr" <<t.sourceAddress << "Dst Addr "<< t.destinationAddress);
NS_LOG_UNCOND("Sent Packets=" <<iter->second.txPackets);
NS_LOG_UNCOND("Received Packets =" <<iter->second.rxPackets);
NS_LOG_UNCOND("Lost Packets =" <<iter->second.txPackets-iter->second.rxPackets);
NS_LOG_UNCOND("Packet delivery ratio =" <<iter->second.rxPackets*100/iter->second.txPackets << "%");
NS_LOG_UNCOND("Packet loss ratio =" << (iter->second.txPackets-iter->second.rxPackets)*100/iter->second.txPackets << "%");
NS_LOG_UNCOND("Delay =" <<iter->second.delaySum);
NS_LOG_UNCOND("Jitter =" <<iter->second.jitterSum);
NS_LOG_UNCOND("Throughput =" <<iter->second.rxBytes * 8.0/(iter->second.timeLastRxPacket.GetSeconds()-iter->second.timeFirstTxPacket.GetSeconds())/1024<<"Kbps");

SentPackets = SentPackets +(iter->second.txPackets);
ReceivedPackets = ReceivedPackets + (iter->second.rxPackets);
LostPackets = LostPackets + (iter->second.txPackets-iter->second.rxPackets);
AvgThroughput = AvgThroughput + (iter->second.rxBytes * 8.0/(iter->second.timeLastRxPacket.GetSeconds()-iter->second.timeFirstTxPacket.GetSeconds())/1024);
Delay = Delay + (iter->second.delaySum);
Jitter = Jitter + (iter->second.jitterSum);

j = j + 1;

}

AvgThroughput = AvgThroughput/j;
NS_LOG_UNCOND("--------Total Results of the simulation----------"<<std::endl);
NS_LOG_UNCOND("Total sent packets  =" << SentPackets);
NS_LOG_UNCOND("Total Received Packets =" << ReceivedPackets);
NS_LOG_UNCOND("Total Lost Packets =" << LostPackets);
NS_LOG_UNCOND("Packet Loss ratio =" << ((LostPackets*100)/SentPackets)<< "%");
NS_LOG_UNCOND("Packet delivery ratio =" << ((ReceivedPackets*100)/SentPackets)<< "%");
NS_LOG_UNCOND("Average Throughput =" << AvgThroughput<< "Kbps");
NS_LOG_UNCOND("End to End Delay =" << Delay);
NS_LOG_UNCOND("End to End Jitter delay =" << Jitter);
NS_LOG_UNCOND("Total Flod id " << j);
monitor->SerializeToXmlFile("manet-routing.xml", true, true);

*/
          //std::cout << "Total Tx bits = " << txByteCounter * 8 << std:: endl;
           //std::cout << "Total Tx packets = " << GetRoutingStats().GetTxPkts() << std:: endl;



         //std::cout << "Total Rx bits = " << bytesTotal* 8 << std:: endl;
        //  std::cout << "Total Rx packets = " << packetsReceived << std:: endl;
  std::cout << "PDR = " << ((double)GetRoutingStats().GetRxPkts())/ ((double)GetRoutingStats().GetTxPkts()) << std:: endl;

  std::cout << "Avrg thput = " << ((GetRoutingStats ().GetRxBytes ())* 8) / (finalSimTime - Seconds(realAppStart)).GetSeconds () / 1000.0 << " kbps" << std:: endl;
  //std::cout << "Avrg thput = " << ReceiveRate << " kbps" << std:: endl;

  Simulator::Destroy ();
//  return 0;
}
