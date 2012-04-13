CREATE TABLE IF NOT EXISTS E_DIFF_CROSS_SECTION (

  name             VARCHAR(30)    NOT NULL,
  measurement_tag  VARCHAR(30)    DEFAULT '0' NOT NULL,

  Sigma            FLOAT          NOT NULL,
  Sigma_units      VARCHAR(40)    DEFAULT '10-38 cm2/GeV/nucleon' NOT NULL,
  dSigma           FLOAT          DEFAULT '0',

  E                FLOAT          NOT NULL,
  E_units          VARCHAR(30)    DEFAULT 'GeV' NOT NULL,

  EP               FLOAT          NOT NULL,
  EP_units         VARCHAR(30)    DEFAULT 'GeV' NOT NULL,

  Theta            FLOAT          NOT NULL,
  Theta_units      VARCHAR(30)    DEFAULT 'degrees' NOT NULL,

  Q2               FLOAT          NOT NULL,
  Q2_units         VARCHAR(30)    DEFAULT 'GeV2' NOT NULL,

  W2               FLOAT          NOT NULL,
  W2_units         VARCHAR(30)    DEFAULT 'GeV2' NOT NULL,

  Nu               FLOAT          NOT NULL,
  Nu_units         VARCHAR(30)    DEFAULT 'GeV' NOT NULL,

  Epsilon          FLOAT          DEFAULT '0',
  Gamma            FLOAT          DEFAULT '0',
  x                FLOAT          DEFAULT '0',

  PRIMARY KEY (name, measurement_tag, E, EP, Theta)
);
