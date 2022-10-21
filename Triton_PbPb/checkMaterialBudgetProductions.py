import ROOT

output_file = ROOT.TFile.Open("MaterialRadiusCheck.root","recreate")
input_eff_r_material = ROOT.TFile.Open("out/EfficiencyHe3_LHC22b9.root")
input_eff_std_material = ROOT.TFile.Open("out/EfficiencyHe3_LHC22b9_2.root")

cent = [[0,5],[5,10],[30,50],[999,999]]
matt = [["A","antihelium"],["M","helium"]]
outlier = [2.15,2.9]

ROOT.gStyle.SetOptStat(0)
gCorrect=[ROOT.TGraphErrors(),ROOT.TGraphErrors()]

for cent_class in cent:
    if cent_class[0]<100:
        for matt_label in matt:
            c = ROOT.TCanvas(f"c{matt_label[0]}Compare_{cent_class[0]}_{cent_class[1]}")
            l = ROOT.TLegend(0.4,0.2,0.85,0.45)
            h_std = input_eff_std_material.Get(f"f{matt_label[0]}Eff_TPC_{cent_class[0]}_{cent_class[1]}")
            h_r = input_eff_r_material.Get(f"f{matt_label[0]}Eff_TPC_{cent_class[0]}_{cent_class[1]}")
            # if cent_class[0]<1:
            #     for bin in outlier:
            #         h_std.SetBinContent(h_std.FindBin(bin),0)
            #         h_std.SetBinError(h_std.FindBin(bin),0)
            #         h_r.SetBinContent(h_r.FindBin(bin),0)
            #         h_r.SetBinError(h_r.FindBin(bin),0)
            h_r.SetLineColor(ROOT.kRed)
            h_r.SetMarkerColor(ROOT.kRed)
            h_r.SetMarkerStyle(20)
            h_r.SetMarkerSize(1.0)
            h_std.SetLineColor(ROOT.kBlack)
            h_std.SetMarkerColor(ROOT.kBlack)
            h_std.SetMarkerStyle(25)
            h_std.SetMarkerSize(1.0)
            h_r.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
            h_std.GetXaxis().SetRangeUser(0.,10.)
            h_r.GetXaxis().SetRangeUser(0.,10.)
            h_std.GetYaxis().SetRangeUser(0.,1.)
            h_r.GetYaxis().SetRangeUser(0.,1.)
            l.AddEntry(h_std,"Standard material budget")
            l.AddEntry(h_r,"R-dependent material budget")
            h_r.Draw()
            h_std.Draw("same")
            l.Draw("same")
            c.Print(f"{c.GetName()}.pdf")

            cRatio = ROOT.TCanvas(f"c{matt_label[0]}Ratio_{cent_class[0]}_{cent_class[1]}")
            h_rr = ROOT.TH1D(h_r)
            h_rr.Divide(h_std)
            h_rr.SetTitle(f"{matt_label[1]} efficiency ratio, {cent_class[0]}-{cent_class[1]}%")
            h_rr.GetYaxis().SetTitle("R-dependent / Standard")
            h_rr.GetYaxis().SetRangeUser(0.5,1.5)
            h_rr.Fit("pol0","R","",2.,10.)
            lH = ROOT.TLine(0.25,1.,10.,1.)
            lH.SetLineColor(ROOT.kBlack)
            lH.SetLineStyle(ROOT.kDashed)
            h_rr.Draw()
            lH.Draw("same")
            cRatio.Print(f"{cRatio.GetName()}.pdf")
            output_file.cd()
            h_rr.Write()
            if matt_label[0]=="A":
                gCorrect[0].AddPoint(0.5*(cent_class[0]+cent_class[1]),h_rr.GetFunction("pol0").GetParameter(0))
                gCorrect[0].SetPointError(gCorrect[0].GetN()-1,0.,h_rr.GetFunction("pol0").GetParError(0))
            else:
                gCorrect[1].AddPoint(0.5*(cent_class[0]+cent_class[1]),h_rr.GetFunction("pol0").GetParameter(0))
                gCorrect[1].SetPointError(gCorrect[1].GetN()-1,0.,h_rr.GetFunction("pol0").GetParError(0))
    else:
        for i_g, graph in enumerate(gCorrect):
            graph.Fit("pol0")
            graph.Write(f"Graph{matt[i_g][0]}")

output_file.Close()