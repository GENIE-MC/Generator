CREATE TABLE IF NOT EXISTS CROSS_SECTION (

  name             VARCHAR(30)   NOT NULL,
  measurement_tag  VARCHAR(30)   DEFAULT '0' NOT NULL,

  xsec             FLOAT         NOT NULL,
  stat_err_p       FLOAT,
  stat_err_m       FLOAT,
  syst_err_p       FLOAT,
  syst_err_m       FLOAT,

  xsec_units       VARCHAR(40)   DEFAULT '10-38 cm2/GeV/nucleon' NOT NULL,
  xsec_norm        VARCHAR(10)   DEFAULT 'E'                     NOT NULL,
  stat_err_type    VARCHAR(10)   DEFAULT 'xsec'                  NOT NULL,
  syst_err_type    VARCHAR(10)   DEFAULT 'xsec'                  NOT NULL,

  E                FLOAT         NOT NULL,
  E_min            FLOAT,
  E_max            FLOAT,
  E_units          VARCHAR(30)   DEFAULT 'GeV' NOT NULL,
  E_frame          VARCHAR(30)   DEFAULT 'lab' NOT NULL,

  PRIMARY KEY (name, measurement_tag, E)
);
