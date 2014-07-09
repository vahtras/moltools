#!/usr/bin/env python
#-*- coding: utf-8 -*-

class Dic:
    def __init__(self):
        rho3 = { 0.00 : [] }
        rho2 = { 0.00 : rho3 }
        rho1 = { 0.00 : rho2 }
        tau = { 0.00 : rho1}
        theta = { 0.00 : tau }
        r = { 0.00 : theta }
        self.dic = r
    def getVal(self, r, theta, tau, rho1, rho2, rho3):
        if self.dic.has_key( r ):
            if self.dic[r].has_key( theta ):
                if self.dic[r][theta].has_key(tau):
                    if self.dic[r][theta][tau].has_key(rho1):
                        if self.dic[r][theta][tau][rho1].has_key(rho2):
                            if self.dic[r][theta][tau][rho1][rho2].has_key(rho3):
                                return self.dic[r][theta][tau][rho1][rho2][rho3]
    def setVal(self, r, theta, tau, rho1, rho2, rho3, val ): #, rho1, rho2, rho3, val):
        tmpr = r
        tmptheta = theta
        tmptau = tau
        tmprho1 = rho1
        tmprho2 = rho2
        tmprho3 = rho3
        if self.dic.has_key( r ):
            if self.dic[r].has_key( theta ):
                if self.dic[r][theta].has_key(tau):
                    if self.dic[r][theta][tau].has_key(rho1):
                        if self.dic[r][theta][tau][rho1].has_key(rho2):
                            if self.dic[r][theta][tau][rho1][rho2].has_key(rho3):
                                self.dic[r][theta][tau][rho1][rho2][rho3] = val
                            else:
                                self.dic[r][theta][tau][rho1][rho2][rho3] = val
                        else:
                            rho3 = { tmprho3 : val }
                            self.dic[r][theta][tau][rho1][rho2] = rho3
                    else:
                        rho3 = { tmprho3 : val }
                        rho2 = { tmprho2 : rho3 }
                        self.dic[r][theta][tau][rho1] = rho2
                else:
                    rho3 = { tmprho3 : val }
                    rho2 = { tmprho2 : rho3 }
                    rho1 = { tmprho1 : rho2 }
                    self.dic[r][theta][tau] = rho1
            else:
                rho3 = { tmprho3 : val }
                rho2 = { tmprho2 : rho3 }
                rho1 = { tmprho1 : rho2 }
                tau = { tmptau : rho1 }
                self.dic[r][theta] = tau
        else:
            rho3 = { tmprho3 : val }
            rho2 = { tmprho2 : rho3 }
            rho1 = { tmprho1 : rho2 }
            tau = { tmptau : rho1 }
            theta = { tmptheta : tau }
            self.dic[r] = theta
