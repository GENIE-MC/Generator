CREATE TABLE IF NOT EXISTS BEAM_FLUX (

  name             VARCHAR(30)   NOT NULL,
  beam             VARCHAR(30)   NOT NULL,

  energy           FLOAT         NOT NULL,
  energy_min       FLOAT         NOT NULL,
  energy_max       FLOAT         NOT NULL,
  energy_units     VARCHAR(10)   DEFAULT 'GeV' NOT NULL,
  energy_frame     VARCHAR(10)   DEFAULT 'lab' NOT NULL,

  flux             FLOAT         NOT NULL,
  dflux_pos        FLOAT         DEFAULT 0 NOT NULL,
  dflux_neg        FLOAT         DEFAULT 0 NOT NULL,
  flux_units       VARCHAR(40)   DEFAULT 'cm-2 sec-1 srad-1 GeV-1' NOT NULL,

  PRIMARY KEY (name, beam, energy)
);
